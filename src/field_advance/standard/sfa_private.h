#ifndef _sfa_private_h_
#define _sfa_private_h_

// Standard field advance implementation

#ifndef IN_sfa
#error "Do not include sfa_private.h; include field_advance.h"
#endif

#define IN_field_advance
#include "../field_advance_private.h"

typedef struct material_coefficient {
  float decayx, drivex;         // Decay of ex and drive of (curl H)x and Jx
  float decayy, drivey;         // Decay of ey and drive of (curl H)y and Jy
  float decayz, drivez;         // Decay of ez and drive of (curl H)z and Jz
  float rmux, rmuy, rmuz;       // Reciprocle of relative permeability
  float nonconductive;          // Divergence cleaning related coefficients
  float epsx, epsy, epsz; 
  float pad[3];                 // For 64-byte alignment and future expansion
} material_coefficient_t;

typedef struct sfa_params {
  material_coefficient_t * mc;
  int n_mc;
  float damp;
} sfa_params_t;

BEGIN_C_DECLS

// In standard_field_advance.c

void
delete_standard_field_array( field_array_t * RESTRICT fa );

void
clear_jf( field_array_t * RESTRICT fa );

void
clear_rhof( field_array_t * RESTRICT fa );

// In advance_b.c

// advance_b applies the following difference equation to the fields:
//   c B_new = c B_old - frac c dt curl E

void
advance_b( field_array_t * RESTRICT fa,
           float                    frac );

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
//
// vacuum_advance_e is the high performance version for uniform regions
//
// FIXME: Currently, frac must be 1.

void
advance_e( field_array_t * RESTRICT fa,
           float                    frac );

void
vacuum_advance_e( field_array_t * RESTRICT fa,
                  float                    frac );

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
//
// vacuum_energy_f is the high performance version for uniform regions

void
energy_f( /**/  double        * RESTRICT en, // 6 elem array
          const field_array_t * RESTRICT fa );

void
vacuum_energy_f( /**/  double        * RESTRICT en, // 6 elem array
                 const field_array_t * RESTRICT fa );

// In compute_curl_b.c

// compute_curl_b applies the following difference equations to the
// fields
//   tca = c dt curl ( c B / mu_r )
// Note: compute_curl_b is structurally the same as advance_e; updates
// to one likely should be replicated in the other.
//
// vacuum_compute_curl_b is the high performance version for uniform regions

void
compute_curl_b( field_array_t * RESTRICT fa );

void
vacuum_compute_curl_b( field_array_t * RESTRICT fa );

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
//
// vacuum_compute_rhob is the high performance version for uniform regions

void
compute_rhob( field_array_t * RESTRICT fa );

void
vacuum_compute_rhob( field_array_t * RESTRICT fa );

// In compute_div_e_err.c

// compute_div_e_err applies the following difference equation:
//   div_e_err = nonconductive [ ( div eps_r E ) - ( rho_f + rho_b )/eps_0 ]
// Note: The structure of this routine is identical to compute_rhob,
// so updates to one should likely be replicated in the other.
//
// vacuum_compute_div_e_err is the high performance version for uniform regions

void
compute_div_e_err( field_array_t * RESTRICT fa );

void
vacuum_compute_div_e_err( field_array_t * RESTRICT fa );

// In compute_rms_div_e_err.c

// compute_rms_div_e_err returns
//   eps0 sqrt( Integral |div_e_err|^2 / Volume )
// This has units of (electric charge / volume).  The integrals are
// done over all the domains.  The volume is the total volume of all
// domains.  Every processor gets the same value.  Note that this
// function does _not_ update or recompute div_e_err.

double
compute_rms_div_e_err( const field_array_t * RESTRICT fa );

// In clean_div_e.c

// clean_div_e applies the following difference equation:
//   E_new = E_old + drive alpha dt grad div_e_err
// div_e_err is not updated or recomputed by this function.
//
// vacuum_clean_div_e is the high performance version for uniform regions

void
clean_div_e( field_array_t * RESTRICT fa );

void
vacuum_clean_div_e( field_array_t * RESTRICT fa );

// In compute_div_b_err.c

// compute_div_b_err applies the following difference equation:
//   div_b_err = div cB

void
compute_div_b_err( field_array_t * RESTRICT fa );

// In compute_rms_div_b_err.c

// compute_rms_div_b_err returns
//   eps0 sqrt( Integral |div_b_err|^2 / Volume )
// This has units of (electric charge / volume) ... yes electric
// charge.  The integrals are done over all the domains. The volume is
// the total volume of all domains.  Every processor gets the same
// value.  Uses the value of div_b_err already stored.  It _does_
// _not_ recompute div_b_err.

double
compute_rms_div_b_err( const field_array_t * RESTRICT fa );

// In clean_div_b.c

// clean_div_b applies the following difference equation:
//   cB_new = cB_old + alpha dt grad div_b_err
// alpha is picked to rapidly reduce the rms_div_b_err

void
clean_div_b( field_array_t * RESTRICT fa );

// Internode functions

// In remote.c

double
synchronize_tang_e_norm_b( field_array_t * RESTRICT fa );

void
synchronize_jf( field_array_t * RESTRICT fa );

void
synchronize_rho( field_array_t * RESTRICT fa );

// In local.c

void
local_ghost_tang_b( field_t      * ALIGNED(128) f,
                    const grid_t *              g );

void
local_ghost_norm_e( field_t      * ALIGNED(128) f,
                    const grid_t *              g );

void
local_ghost_div_b( field_t      * ALIGNED(128) f,
                   const grid_t *              g );

void
local_adjust_tang_e( field_t      * ALIGNED(128) f,
                     const grid_t *              g );

void
local_adjust_div_e( field_t      * ALIGNED(128) f,
                    const grid_t *              g );

void
local_adjust_norm_b( field_t      * ALIGNED(128) f,
                     const grid_t *              g );

void
local_adjust_jf( field_t      * ALIGNED(128) f,
                 const grid_t *              g );

void
local_adjust_rhof( field_t      * ALIGNED(128) f,
                   const grid_t *              g );

void
local_adjust_rhob( field_t      * ALIGNED(128) f,
                   const grid_t *              g );

// In remote.c

void
begin_remote_ghost_tang_b( field_t      * ALIGNED(128) f,
                           const grid_t *              g );

void
end_remote_ghost_tang_b( field_t      * ALIGNED(128) f,
                         const grid_t *              g );

void
begin_remote_ghost_norm_e( field_t      * ALIGNED(128) f,
                           const grid_t *              g );

void
end_remote_ghost_norm_e( field_t      * ALIGNED(128) f,
                         const grid_t *              g );

void
begin_remote_ghost_div_b( field_t      * ALIGNED(128) f,
                          const grid_t *              g );

void
end_remote_ghost_div_b( field_t      * ALIGNED(128) f,
                        const grid_t *              g );

END_C_DECLS

#endif // _sfa_private_h_
