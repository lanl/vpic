#ifndef _field_h_
#define _field_h_

// FIXME: THE HOST PROCESSED FIELD KERNELS SHOULD BE UPDATED TO USE
// SCALAR FMA INSTRUCTIONS WITH COMMENSURATE ROUND-OFF PROPERTIES TO
// THE FMA INSTRUCTIONS USED ON THE PIPELINE PROCESSED FIELDS!

// FIXME: (nx>1) ? (1/dx) : 0 TYPE LOGIC SHOULD BE FIXED SO THAT NX
// REFERS TO THE GLOBAL NUMBER OF CELLS IN THE X-DIRECTION (NOT THE
// _LOCAL_ NUMBER OF CELLS).  THIS LATENT BUG IS NOT EXPECTED TO
// AFFECT ANY PRACTICAL SIMULATIONS.

// FIXME: EXTERNAL DIAGNOSTICS THAT READ THESE WILL NEED TO BE UPDATED
// TO REFLECT SPU USAGE ALIGNMENT CHANGES.
  
#include <material.h>

// This module implements the following the difference equations on a
// superhexahedral domain decomposed Yee-mesh:
//  
// advance_b -> Finite Differenced Faraday
//   cB_new = cB_old - frac c dt curl E
//
// advance_e -> Exponentially Differenced Ampere
//   TCA_new = -damp TCA_old + (1+damp) c dt curl (cB/mu_r)
//   E_new = decay E_old + drive [ TCA_new - (dt/eps0) J  ]
//     where the coefficients decay and drive are defined by the
//     diagonal materials properties tensor. That is:
//       decay_xx = exp(-alpha_xx)
//       drive_xx = [ 1 - exp(-alpha_xx) ] / (alpha_xx epsr_xx)
//       alpha_xx = sigma_xx dt / (eps0 epsr_xx)
//
// div_clean_b -> Marder pass on magnetic fields
//   cB_new = cB_old + D dt grad div cB_old
//     where the diffusion coefficient D is automatically selected to
//     rapidly reduce RMS divergence error assuming divergences errors
//     are due to accumulation of numerical roundoff when integrating
//     Faraday. See clean_div.c for details.
//     
// div_clean_e -> Modified Marder pass on electric fields
//   E_new = E_old + drive D dt grad err_mul div ( epsr E_old - rho/eps0 )
//     Since the total rho may not be known everywhere (for example in
//     conductive regions), err_mul sets the divergence error to zero
//     in regions where total rho is unknown. Same considerations as
//     div_clean_b for the diffusion coefficient determination. See
//     clean_div.c for details.
//
// Domain fields are stored in a field_t array which is FORTRAN
// indexed 0:nx+1,0:ny+1,0:nz+1.  The fields are located on a
// staggered Yee-mesh
//
//   f(i,j,k).ex  @ i+0.5,j,k (all 1:nx,1:ny+1,1:nz+1; int 1:nx,2:ny,2:nz)
//   f(i,j,k).ey  @ i,j+0.5,k (all 1:nx+1,1:ny,1:nz+1; int 2:nx,1:ny,2:nz)
//   f(i,j,k).ez  @ i,j,k+0.5 (all 1:nx+1,1:ny+1,1:nz; int 2:nx,2:ny,1:nz)
//   f(i,j,k).cbx @ i,j+0.5,k+0.5 (all 1:nx+1,1:ny,1:nz; int 2:nx,1:ny,1:nz)
//   f(i,j,k).cby @ i+0.5,j,k+0.5 (all 1:nx,1:ny+1,1:nz; int 1:nx,2:ny,1:nz)
//   f(i,j,k).cbz @ i+0.5,j+0.5,k (all 1:nx,1:ny,1:nz+1; int 1:nx,1:ny,2:nz)
//   f(i,j,k).rhof @ i,j,k (all 1:nx+1,1:ny+1,1:nz+1; int 2:nx,2:ny,2:nz)
//   f(i,j,k).div_b_err(i,j,k) @ i+0.5,j+0.5,k+0.5 (all+int 1:nx,1:ny,1:nz)
//   f(i,j,k).rhob, div_e_err, nmat is on same mesh as rhof
//   f(i,j,k).jfx, jfy, jfz is on same mesh as ex,ey,ez respectively
//   f(i,j,k).tcax ,tcay, tcaz is on same mesh ex,ey,ez respectively
//   f(i,j,k).ematx, ematy, ematz are on same mesh as ex,ey,ez respectively
//   f(i,j,k).fmatx, fmaty, fmatz are on same mesh as cbx,cby,cbz respectively
//   f(i,j,k).cmat is on same mesh as div_b_err
//
// Alternatively, ex,ey,ez / jfx,jfy,jfz / tcax,tcay,tcaz /
// ematx,ematy,ematz are all on the "edge mesh". cbx,cby,cbz /
// fmatx,fmaty,fmatz are all on the "face
// mesh". rhof,rhob,div_e_err,nmat are on the "nodes mesh".
// div_b_err,cmat are on the "cell mesh".
// 
// Above, for "edge mesh" quantities, interior means that the
// component is not a tangential field directly on the surface of the
// domain. For "face mesh" quantities, interior means that the
// component is not a normal field directly on the surface of the
// computational domain. For "node mesh" quantities, interior means
// that the component is not on the surface of the computational
// domain. All "cell mesh" quantities are interior to the computation
// domains.
//
// Boundary conditions and field communications are implemented in
// part by ghost values:
//   advance_e uses tangential b ghosts
//     cbx => 1:nx+1,{0|ny+1},1:nz and 1:nx+1,1:ny,{0|nz+1}
//     cby => 1:nx,1:ny+1,{0|nz+1} and {0|nx+1},1:ny+1,1:nz
//     cbz => {0|nx+1},1:ny,1:nz+1 and 1:nx,{0|ny+1},1:nz+1
//   advance_b does not use any ghosts
//   clean_div_e uses normal e ghosts.
//   clean_div_b uses derr ghosts.
//
// To setup the field solver, the following sequence is suggested:
//   grid_t grid[1];
//   material_t *material_list = NULL;
//   fields_t * ALIGNED(16) fields = NULL;
//   material_coefficient_t * ALIGNED(16) material_coefficients = NULL;
//
//   grid->nx = 32; grid->ny = 32; grid->nz = 32;
//   ...
//   new_material("Vacuum",1,1,0,&material_list);
//   ...
//   material_coefficients = new_material_coefficients(grid,material_list);
//   fields = new_fields(grid);
// 
//   ... Set the initial field values and place materials ...
//
//   synchronize_fields(fields,grid);
//
// Note: synchronize_fields makes sure shared faces on each domain
// have the same fields (in case the user messed up setting the
// initial fields or errors in the source terms or different floating
// point properties on different nodes cause the shared faces to have
// different fields).
//   
// To advance the fields in a PIC simulation with TCA radation damping
// and periodic divergence cleaning, the following sequence is
// suggested:
//   ... push/accumulate particles using E_0, B_0 => Jf_0.5
//   advance_b( fields, grid, 0.5  );                  => B_0 to B_0.5
//   advance_e( fields, material_coefficients, grid ); => E_0 to E_1
//   advance_b( fields, grid, 0.5  );                  => B_0.5 to B_1
//   if( should_clean_div_e ) {
//     ... adjust rho_f, rho_b and/or rho_c as necessary
//     do {
//       rms_err = clean_div_e( fields, material_coefficients, grid ); 
//     } while( rms_err_too_high );
//   }
//   if( should_clean_div_b ) {
//     do {
//       rms_err = clean_div_b( fields, grid );
//     } while( rms_err_too_high );
//   }
//   if( should_sync_fields ) synchronize_fields(fields,grid);
//
// To clean up the field solver, the following sequence is suggested:
//   delete_materials(&materials_list);
//   delete_materials_coefficients(&materials_coefficients);
//   delete_fields(&fields);

// Field arrays shall be a (nx+2) x (ny+2) x (nz+2) allocation indexed
// FORTRAN style from (0:nx+1,0:ny+1,0:nz+1).  Fields for voxels on
// the surface of the local domain (for example h(0,:,:) or
// h(nx+1,:,:)) are used to store values of the field array on
// neighbor processors or used to enforce bonundary conditions.

typedef struct field {
  float ex, ey, ez;                // Electric field
  float cbx, cby, cbz;             // Magnetic field
  float tcax, tcay, tcaz;          // TCA fields
  float jfx, jfy, jfz;             // Free current
  material_id ematx, ematy, ematz; // Material at edge centers
  material_id fmatx, fmaty, fmatz; // Material at face centers
  material_id nmat, cmat;          // Material at cell centers and nodes
  float rhof, rhob;                // Free and bound charge density
  float div_e_err, div_b_err;      // Divergence errors
# ifdef USE_CELL_SPUS
  float _pad[12];                  // 128-byte align (next power of two)
# else
  /**/                             // 16-byte align
# endif
} field_t;

// Hydro arrays shall be a (nx+2) x (ny+2) x (nz+2) allocation indexed
// FORTRAN style from (0:nx+1,0:ny+1,0:nz+1).  Hydros for voxels on
// the surface of the local domain (for example h(0,:,:) or
// h(nx+1,:,:)) are not used.

typedef struct hydro {
  float rho;           // Charge density         => < q f >
  float jx, jy, jz;    // Current density        => < q v_i f >
  float ke;            // Kinetic energy density => < m c^2 (gamma-1) f >
  float px, py, pz;    // Momentum density       => < p_i f >
  float txx, tyy, tzz; // Stress diagonal        => < p_i v_j f >, i==j
  float tyz, tzx, txy; // Stress off-diagonal    => < p_i v_j f >, i!=j
# ifdef USE_CELL_SPUS
  float _pad[2];       // 64-byte align (next power of two)
# else
  float _pad[2];       // 16-byte align
# endif
} hydro_t;

// Interpolator arrays shall be a (nx+2) x (ny+2) x (nz+2) allocation
// indexed FORTRAN style from (0:nx+1,0:ny+1,0:nz+1). Interpolators
// for voxels on the surface of the local domain (for example
// fi(0,:,:) or fi(nx+1,:,:)) are not used.

typedef struct interpolator {
  float ex, dexdy, dexdz, d2exdydz;
  float ey, deydz, deydx, d2eydzdx;
  float ez, dezdx, dezdy, d2ezdxdy;
  float cbx, dcbxdx;
  float cby, dcbydy;
  float cbz, dcbzdz;
# ifdef USE_CELL_SPUS
  float _pad[14]; // 128-byte align (next power of two)
# else
  float _pad[2];  // 16-byte align
# endif
} interpolator_t;

// Accumulator arrays shall be a (nx+2)x(ny+2)x(nz+2)x(1+n_pipeline)
// allocation indexed FORTRAN style.  That is, the accumulator array
// is a 4d array.  a(:,:,:,0) is the accumulator used by the host
// processor.  a(:,:,:,1:n_pipeline) are the accumulators used by
// pipelines during operations.  Like the interpolator, accumualtors
// on the surface of the local domain are not used.

typedef struct accumulator {
  float jx[4];   // jx0@(0,-1,-1),jx1@(0,1,-1),jx2@(0,-1,1),jx3@(0,1,1)
  float jy[4];   // jy0@(-1,0,-1),jy1@(-1,0,1),jy2@(1,0,-1),jy3@(1,0,1)
  float jz[4];   // jz0@(-1,-1,0),jz1@(1,-1,0),jz2@(-1,1,0),jz3@(1,1,0)
# ifdef USE_CELL_SPUS
  float _pad[4]; // 64-byte align (next power of two )
# else
  /**/           // 16-byte align
# endif
} accumulator_t;

BEGIN_C_DECLS

// Constructors and destructors

// In field_structors.c

field_t * ALIGNED(16)
new_field( const grid_t * g );

void
delete_field( field_t ** ALIGNED(16) f );

void
clear_field( field_t      * ALIGNED(16) f,
             const grid_t *             g );

void
clear_jf( field_t      * ALIGNED(16) f,
          const grid_t *             g );

void
clear_rhof( field_t      * ALIGNED(16) f,
            const grid_t *             g );

hydro_t * ALIGNED(16)
new_hydro( const grid_t * g );

void
delete_hydro( hydro_t ** ALIGNED(16) h );

void
clear_hydro( hydro_t      * ALIGNED(16) h,
             const grid_t *             g );

interpolator_t * ALIGNED(128)
new_interpolator( const grid_t * g );

void
delete_interpolator( interpolator_t ** ALIGNED(128) fi );

void
clear_interpolator( interpolator_t * ALIGNED(128) fi,
                    const grid_t   *              g );

accumulator_t * ALIGNED(128)
new_accumulators( const grid_t * g );

void
delete_accumulators( accumulator_t ** ALIGNED(128) a );

// Time stepping functions

// In load_interpolator.c

// Going into load_interpolator, the field array f contains the
// current information such that the fields can be interpolated to
// particles within the local domain.  Load interpolate computes the
// field array into a set of interpolation coefficients for each voxel
// inside the local domain suitable for use by the particle update
// functions.

void
load_interpolator( interpolator_t * ALIGNED(128) fi,
                   const field_t  * ALIGNED(16)  f,
                   const grid_t   *              g );


// In clear_accumulators.c

// This zeros out all the accumulator arrays in a pipelined fashion.

void
clear_accumulators( accumulator_t * ALIGNED(128) a,
                    const grid_t  *              g );

// In reduce_accumulators.c

// Going into reduce_accumulators, the host processor and the pipeline
// processors have each accumulated values to their personal
// accumulators.  This reduces the pipeline accumulators into the host
// accumulator with a pipelined horizontal reduction (a deterministic
// reduction).

void
reduce_accumulators( accumulator_t * ALIGNED(128) a,
                     const grid_t  *              g );

// In unload_accumulator.c

// Going into unload_accumulator, the accumulator contains 4 times the
// net amount of charge that crossed the quarter face associated with
// each accumulator component (has units of physical charge, i.e. C)
// computed by the advance_p functions.  unload_accumulator computes
// the physical current density (A/m^2 in MKS units) associated with
// all local quarter faces and accumulates the local quarter faces to
// local field array jf.  unload_accumulator assumes all the pipeline
// accumulators have been reduced into the host accumulator.

void
unload_accumulator( field_t             * ALIGNED(16)  f, 
                    const accumulator_t * ALIGNED(128) a,
                    const grid_t        *              g );

// In advance_b.c

// advance_b applies the following difference equation to the fields:
//   c B_new = c B_old - frac c dt curl E

void
advance_b( field_t      * ALIGNED(16) f,
           const grid_t *             g,
           float                      frac );

// In advance_e.c

// advance_e applies the following difference equations to the fields
//   tca_new = ( 1 + damp ) c dt curl ( c B / mu_r ) -
//               damp tca_old
//   E_new = decay E_old + drive [ tca_new - (dt/eps0) Jf ]
// where: 
//   damp is numerical Cherenkov damping parameter
//   decay = exp( -alpha )
//   drive = ( 1 - decay ) / ( alpha eps_r )
//     alpha = sigma dt / ( eps_0 eps_r )
//     sigma is the conductivity (J_c = sigma E)
// Note: advance_e is structurally the same as compute_curl_b.
// Updates to one likely should be replicated in the other.

void
advance_e( field_t                      * ALIGNED(16) f,
           const material_coefficient_t * ALIGNED(16) m,
           const grid_t                 *             g );

// Diagnostic functions

// In energy_f.c

// This computes 6 components of field energy of the system.  The
// return value is:
//   energy[0] = integral_volume 0.5 Dx Ex d^3 r ... Ex field energy
//   energy[1] = integral_volume 0.5 Dy Ey d^3 r ... Ey field energy
//   energy[2] = integral_volume 0.5 Dz Ez d^3 r ... Ez field energy
//   energy[3] = integral_volume 0.5 Hx Bx d^3 r ... Bx field energy
//   energy[4] = integral_volume 0.5 Hy By d^3 r ... By field energy
//   energy[5] = integral_volume 0.5 Hz Bz d^3 r ... Bz field energy
// Thus, sum energy = 0.5 integral_volume ( D.E + H.B ) d^3 r
// All nodes get the same result.

// FIXME: Should this function perform the global reduce?

void
energy_f( double                       *             energy, // 6 elem array
          const field_t                * ALIGNED(16) f,
          const material_coefficient_t * ALIGNED(16) m,
          const grid_t                 *             g );

// In compute_curl_b.c

// compute_curl_b applies the following difference equations to the
// fields
//   tca = c dt curl ( c B / mu_r )
// Note: compute_curl_b is structurally the same as advance_e; updates
// to one likely should be replicated in the other.

void
compute_curl_b( field_t                      * ALIGNED(16) f,
                const material_coefficient_t * ALIGNED(16) m,
                const grid_t                 *             g );

// Divergence cleaning functions

// The theory behind the Marder correction is that the Ampere and
// Faraday equations can be modified as follows:
//   pB/pt = -curl E    --> pB/pt = -curl E     + alpha grad div B
//   pD/pt = curl H - J --> pD/pt =  curl H - J + alpha grad ( div D - rho )
// Taking the divergence of the modified equations yield:
//   p(div B)/pt     = alpha laplacian div B
//   p(div D-rho)/pt = alpha laplacian ( div D - rho )
// Since these are sourceless diffusion equation, asymptotically,
//   div B       --> 0 
//   div D - rho --> 0
// In particular, Fourier transforming div B in space shows that a
// given mode decays as exp(-alpha k^2 t). The diffusion coefficient
// alpha controls how fast the divergence goes to zero. Note the long
// wavelength modes decay more slowly and k=0 does not decay at
// all. This is not a problem though owing to properties of the
// divergence operator:
//   div B @ k=0       -> integral over all space of div B -> 0
//   div D - rho @ k=0 -> integral over all space of div D - rho
//                     -> -net charge in universe -> 0
//
// Since div B and div D-rho is ideally zero, the modified equations
// do not change _any_ physics. Further, if for any reason a non-zero
// div B or (div D - rho) occurs, the above modification will drive
// the error back to zero.
//   
// To understand how use this in a simulation, consider the standard
// field update equations for Bx on a Yee mesh without the additional
// term:
//   cBx(1/2) = cBx(-1/2) - (c*dt)(curl E)_x
// Because of finite precision arithmetic, cBx(1/2) can be off with a
// relative error on order the machine's floating point precision, eps
// (~1.2e-7 for IEEE single precision). Over many time steps, these
// errors accumulate. The accumulation process can be thought of as a
// random walk with a RMS step size on the order of ~0.5 eps
// |cBx|. Thus, over a large number of time steps Nt, the PDF of the
// error in cBx for an arbitrary grid point will be closely
// approximated by a Gaussian with zero mean and standard deviation
// ~0.5 eps |cBx| sqrt(Nt). The same holds true for cBy and cBz.
// 
// If it is assumed that the errors between different grid points are
// uncorrelated (a _very_ accurate assumption except for very
// specially prepared field configurations), then the power in various
// spectral modes of div B on a periodic mesh can be shown to be:
//   |div cB(kx,ky,kz)_unclean|^2 ~
//     [eps^2 |cB|^2 Nt/(Nx Ny Nz)][ (sin(pi kx/Nx)/dx)^2 +
//                                   (sin(pi ky/Ny)/dy)^2 +
//                                   (sin(pi kz/Nz)/dz)^2 ]
// To reduce this error accumulation, the grad div B term is applied
// using forward differencing in time (this is the usual Marder pass
// ... strictly local operations, easy and efficient to implement in
// parallel):
//   cBx(1/2)_clean = cBx(1/2)_unclean + 
//       alpha dt grad div cBx(1/2)_unclean
// The power in various modes of cBx(1/2)_clean can be shown to be:
//  |div cB(kx,ky,kz)_clean|^2 ~
//     |div cB(kx,ky,kz)_unclean|^2 
//       { 1 - (4*alpha*dt/dg^2) [ (dg sin(pi kx/Nx)/dx)^2 +
//                                 (dg sin(pi ky/Ny)/dy)^2 +
//                                 (dg sin(pi kz/Nz)/dz)^2 ] }^2
// where dg^-2 = dx^-2 + dy^-2 + dz^-2.
//
// Inspecting the above, if 0 <= alpha dt < dg^2/2, then no component
// of div cB(kx,ky,kz) grows and the divergence cleaning pass is
// numerically stable. Note: This is the same stability criterion as
// the forward differenced diffusion equation.
// 
// If alpha dt = dg^2/4, then shortest wavelength component of div cB
// will be zeroed. Since this is where most of the divergence errors
// are located, this is a relatively good choice.
// 
// If we want to minimize the total RMS divergence error, it can be
// shown (using Parseval's theorem) that the best choice of alpha dt
// on large cubic periodic meshes is:
//   alpha dt ~ 0.388888889 dg^2
// This value is pretty close to optimal on other meshes also. Using
// this value will take the total RMS divergence error to ~0.304 of
// the original value.
// 
// If we assume future contributions to the divergence error are
// uncorrelated with previous contributions (a very accurate
// assumption) and we are only going to clean every Nc time steps,
// then the maximum relative RMS divergence error, div_max will obey
// the following relation asymptotically:
//   div_max^2 = (0.304 div_max)^2    +    0.25 eps^2 Nc
//                      |                        |
//                      |                        |
//       Error left over from previous clean     |
//                        Error accumulated since previous clean
// Solving for Nc yields:
//   Nc ~ 3.63 (div_max/eps)^2
//
// Example:
//   For div_max ~ 1e-6 in single precision, a divergence clean
//   for cB should be done every 255 time steps.
//
// For the clean_div_e, there are two additional considerations. Since
// advance_e uses exponential differencing in time, the diffusion pass
// should be modified to be consistent with advance_e. If divergence
// cleaning were done as part of the time step:
//   E_new = decay E_old + drive [TCA_new - (dt/eps0) J +
//                                (dt/eps0)alpha grad(div eps0 epsr E - rho)]
// Extracting out the divergence cleaning correction yields:
//   E_clean = E_unclean + drive alpha dt grad (div epsr E - rho/eps0)
// Second, in conductive medium, the total charge is not known
// (conduction charge density is not directly computed). Thus, in
// these regions, the diveregence error is considered to be zero. This
// gives the complete modified Marder pass:
//  E_clean = E_unclean +
//            drive alpha dt grad nonconductive (div epsr E - rho/eps0)
 
// In compute_rhob.c

// compute_rhob applies the following difference equation:
//   rho_b = eps0 nonconductive ( div epsr E - rho_f/eps0 )
// rho_b is not computed correctly on absorbing boundary surfaces but
// it does not matter as divergence cleaning is not applied there.
// Note: The structure of this routine is identical to
// compute_div_e_err, so updates to one should likely be replicated in
// the other.

void
compute_rhob( field_t                      * ALIGNED(16) f,
              const material_coefficient_t * ALIGNED(16) m,
              const grid_t                 *             g ); 

// In compute_div_e_err.c

// compute_div_e_err applies the following difference equation:
//   div_e_err = nonconductive [ ( div eps_r E ) - ( rho_f + rho_b )/eps_0 ]
// Note: The structure of this routine is identical to compute_rhob,
// so updates to one should likely be replicated in the other.

void
compute_div_e_err( field_t                      * ALIGNED(16) f,
                   const material_coefficient_t * ALIGNED(16) m,
                   const grid_t                 *             g );

// In compute_rms_div_e_err.c

// compute_rms_div_e_err returns
//   eps0 sqrt( Integral |div_e_err|^2 / Volume )
// This has units of (electric charge / volume).  The integrals are
// done over all the domains.  The volume is the total volume of all
// domains.  Every processor gets the same value.  Note that this
// function does _not_ update or recompute div_e_err.

// FIXME: Should this function do the global reduce?

double
compute_rms_div_e_err( field_t      * ALIGNED(16) f,
                       const grid_t *             g );

// In clean_div_e.c

// clean_div_e applies the following difference equation:
//   E_new = E_old + drive alpha dt grad div_e_err
// div_e_err is not updated or recomputed by this function.

void
clean_div_e( field_t                      * ALIGNED(16) f,
             const material_coefficient_t * ALIGNED(16) m,
	     const grid_t                 *             g );

// In compute_div_b_err.c

// compute_div_b_err applies the following difference equation:
//   div_b_err = div cB

void
compute_div_b_err( field_t      * ALIGNED(16) f,
                   const grid_t *             g );

// In compute_rms_div_b_err.c

// compute_rms_div_b_err returns
//   eps0 sqrt( Integral |div_b_err|^2 / Volume )
// This has units of (electric charge / volume) ... yes electric
// charge.  The integrals are done over all the domains. The volume is
// the total volume of all domains.  Every processor gets the same
// value.  Uses the value of div_b_err already stored.  It _does_
// _not_ recompute div_b_err.

// FIXME: Should this function do the global reduction?

double
compute_rms_div_b_err( field_t      * ALIGNED(16) f,
                       const grid_t *             g );

// In clean_div_b.c

// clean_div_b applies the following difference equation:
//   cB_new = cB_old + alpha dt grad div_b_err
// alpha is picked to rapidly reduce the rms_div_b_err

void
clean_div_b( field_t      * ALIGNED(16) f,
             const grid_t *             g );

// Internode functions

// In remote.c

double
synchronize_tang_e_norm_b( field_t      * ALIGNED(16) f,
                           const grid_t *             g );

void
synchronize_jf( field_t      * ALIGNED(16) f,
                const grid_t *             g );

void
synchronize_rhof( field_t      * ALIGNED(16) f,
                  const grid_t *             g );

void
synchronize_rhob( field_t      * ALIGNED(16) f,
                  const grid_t *             g );

void
synchronize_hydro( hydro_t      * ALIGNED(16) hydro,
                   const grid_t *             g );
             
// INTERNAL USE ONLY functions

// In distribute_voxels.c

// Given a block of voxels to be processed, determine the number of
// voxels and the first voxel a particular job assigned to a pipeline
// should process.  The return voxel is the number of voxels to
// process.
//
// It is assumed that the pipelines will process voxels in FORTRAN
// ordering (e.g. inner loop increments x-index).
//
// jobs are indexed from 0 to n_job-1.  jobs _always_ assign voxels in
// complete bundles of 4 voxels to facillitate vector processing.  If
// job is set to n_job, this function will compute the voxel index of
// the first voxel in the final incomplete bundle and return the
// number of voxels in the final incomplete bundle.

int
distribute_voxels( int x0,  int x1,    // range of x-indices (inclusive)
                   int y0,  int y1,    // range of y-indices (inclusive)
                   int z0,  int z1,    // range of z-indices (inclusive)
                   int job, int n_job, // job ... on [0,n_job-1]
                   int * _x, int * _y, int * _z ); // first voxel to process

// In local.c

void
local_ghost_tang_b( field_t      * ALIGNED(16) f,
                    const grid_t *             g );

void
local_ghost_norm_e( field_t      * ALIGNED(16) f,
                    const grid_t *             g );

void
local_ghost_div_b( field_t      * ALIGNED(16) f,
                   const grid_t *             g );

void
local_adjust_tang_e( field_t      * ALIGNED(16) f,
                     const grid_t *             g );

void
local_adjust_div_e( field_t      * ALIGNED(16) f,
                    const grid_t *             g );

void
local_adjust_norm_b( field_t      * ALIGNED(16) f,
                     const grid_t *             g );

void
local_adjust_jf( field_t      * ALIGNED(16) f,
                 const grid_t *             g );

void
local_adjust_rhof( field_t      * ALIGNED(16) f,
                   const grid_t *             g );

void
local_adjust_rhob( field_t      * ALIGNED(16) f,
                   const grid_t *             g );

void
local_adjust_hydro( hydro_t      * ALIGNED(16) h,
                    const grid_t *             g );

// In remote.c

void
begin_remote_ghost_tang_b( field_t      * ALIGNED(16) f,
                           const grid_t *             g );

void
end_remote_ghost_tang_b( field_t      * ALIGNED(16) f,
                         const grid_t *             g );

void
begin_remote_ghost_norm_e( field_t      * ALIGNED(16) f,
                           const grid_t *             g );

void
end_remote_ghost_norm_e( field_t      * ALIGNED(16) f,
                         const grid_t *             g );

void
begin_remote_ghost_div_b( field_t      * ALIGNED(16) f,
                          const grid_t *             g );

void
end_remote_ghost_div_b( field_t      * ALIGNED(16) f,
                        const grid_t *             g );

END_C_DECLS

#endif // _field_h_
