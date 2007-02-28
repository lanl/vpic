/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version (some data structures based on earlier
 *                    V4PIC versions) 
 *
 */

#ifndef _field_h_
#define _field_h_

#include <common.h>
#include <grid.h>     /* For grid_t */
#include <material.h> /* For material_coefficient_t */

/*****************************************************************************
 * This module implements the following the difference equations on a
 * superhexahedral domain decomposed Yee-mesh:
 *
 * advance_b -> Finite Differenced Faraday
 *   cB_new = cB_old - frac c dt curl E
 *
 * advance_e -> Exponentially Differenced Ampere
 *   TCA_new = -damp TCA_old + (1+damp) c dt curl (cB/mu_r)
 *   E_new = decay E_old + drive [ TCA_new - (dt/eps0) J  ]
 *     where the coefficients decay and drive are defined by the diagonal
 *     materials properties tensor. That is:
 *       decay_xx = exp(-alpha_xx)
 *       drive_xx = [ 1 - exp(-alpha_xx) ] / (alpha_xx epsr_xx)
 *       alpha_xx = sigma_xx dt / (eps0 epsr_xx)
 *
 * div_clean_b -> Marder pass on magnetic fields
 *   cB_new = cB_old + D dt grad div cB_old
 *     where the diffusion coefficient D is automatically selected to rapidly
 *     reduce RMS divergence error assuming divergences errors are due to
 *     accumulation of numerical roundoff when integrating Faraday. See
 *     clean_div.c for details.
 *     
 * div_clean_e -> Modified Marder pass on electric fields
 *   E_new = E_old + drive D dt grad err_mul div ( epsr E_old - rho/eps0 )
 *     Since the total rho may not be known everywhere (for example in 
 *     conductive regions), err_mul sets the divergence error to zero in
 *     regions where total rho is unknown. Same considerations as div_clean_b
 *     for the diffusion coefficient determination. See clean_div.c for
 *     details.
 *
 * Domain fields are stored in a field_t array which is FORTRAN indexed
 *   0:nx+1,0:ny+1,0:nz+1
 * The fields are located on a staggered Yee-mesh
 *   f(i,j,k).ex  @ i+0.5,j,k (all 1:nx,1:ny+1,1:nz+1; int 1:nx,2:ny,2:nz)
 *   f(i,j,k).ey  @ i,j+0.5,k (all 1:nx+1,1:ny,1:nz+1; int 2:nx,1:ny,2:nz)
 *   f(i,j,k).ez  @ i,j,k+0.5 (all 1:nx+1,1:ny+1,1:nz; int 2:nx,2:ny,1:nz)
 *   f(i,j,k).cbx @ i,j+0.5,k+0.5 (all 1:nx+1,1:ny,1:nz; int 2:nx,1:ny,1:nz)
 *   f(i,j,k).cby @ i+0.5,j,k+0.5 (all 1:nx,1:ny+1,1:nz; int 1:nx,2:ny,1:nz)
 *   f(i,j,k).cbz @ i+0.5,j+0.5,k (all 1:nx,1:ny,1:nz+1; int 1:nx,1:ny,2:nz)
 *   f(i,j,k).rhof @ i,j,k (all 1:nx+1,1:ny+1,1:nz+1; int 2:nx,2:ny,2:nz)
 *   f(i,j,k).div_b_err(i,j,k) @ i+0.5,j+0.5,k+0.5 (all+int 1:nx,1:ny,1:nz)
 *   f(i,j,k).rhob, div_e_err, nmat is on same mesh as rhof
 *   f(i,j,k).jfx, jfy, jfz is on same mesh as ex,ey,ez respectively
 *   f(i,j,k).tcax ,tcay, tcaz is on same mesh ex,ey,ez respectively
 *   f(i,j,k).ematx, ematy, ematz are on same mesh as ex,ey,ez respectively
 *   f(i,j,k).fmatx, fmaty, fmatz are on same mesh as cbx,cby,cbz respectively
 *   f(i,j,k).cmat is on same mesh as div_b_err
 *
 * Alternatively, ex,ey,ez / jfx,jfy,jfz / tcax,tcay,tcaz / ematx,ematy,ematz
 * are all on the "edge mesh". cbx,cby,cbz / fmatx,fmaty,fmatz are all on the
 * "face mesh". rhof,rhob,div_e_err,nmat are on the "nodes mesh".
 * div_b_err,cmat are on the "cell mesh".
 * 
 * Above, for "edge mesh" quantities, interior means that the component is not
 * a tangential field directly on the surface of the domain. For "face mesh"
 * quantities, interior means that the component is not a normal field
 * directly on the surface of the computational domain. For "node mesh"
 * quantities, interior means that the component is not on the surface of the
 * computational domain. All "cell mesh" quantities are interior to the
 * computation domains.
 *
 * Boundary conditions and field communications are implemented in part by
 * ghost values:
 *   advance_e uses tangential b ghosts
 *     cbx => 1:nx+1,{0|ny+1},1:nz and 1:nx+1,1:ny,{0|nz+1}
 *     cby => 1:nx,1:ny+1,{0|nz+1} and {0|nx+1},1:ny+1,1:nz
 *     cbz => {0|nx+1},1:ny,1:nz+1 and 1:nx,{0|ny+1},1:nz+1
 *   advance_b does not use any ghosts
 *   clean_div_e uses normal e ghosts.
 *   clean_div_b uses derr ghosts.
 *
 * To setup the field solver, the following sequence is suggested:
 *   grid_t grid[1];
 *   material_t *material_list = NULL;
 *   fields_t * ALIGNED fields = NULL;
 *   material_coefficient_t * ALIGNED material_coefficients = NULL;
 *
 *   grid->nx = 32; grid->ny = 32; grid->nz = 32;
 *   ...
 *   new_material("Vacuum",1,1,0,&material_list);
 *   ...
 *   material_coefficients = new_material_coefficients(grid,material_list);
 *   fields = new_fields(grid);
 * 
 *   ... Set the initial field values and place materials ...
 *
 *   synchronize_fields(fields,grid);
 *
 * Note: synchronize_fields makes sure shared faces on each domain have the
 * same fields (in case the user messed up setting the initial fields or 
 * errors in the source terms or different floating point properties on 
 * different nodes cause the shared faces to have different fields).
 *   
 * To advance the fields in a PIC simulation with TCA radation damping and
 * periodic divergence cleaning, the following sequence is suggested:
 *   ... push/accumulate particles using E_0, B_0 => Jf_0.5
 *   advance_b( fields, grid, 0.5  );                  => B_0 to B_0.5
 *   advance_e( fields, material_coefficients, grid ); => E_0 to E_1
 *   advance_b( fields, grid, 0.5  );                  => B_0.5 to B_1
 *   if( should_clean_div_e ) {
 *     ... adjust rho_f, rho_b and/or rho_c as necessary
 *     do {
 *       rms_err = clean_div_e( fields, material_coefficients, grid ); 
 *     } while( rms_err_too_high );
 *   }
 *   if( should_clean_div_b ) {
 *     do {
 *       rms_err = clean_div_b( fields, grid );
 *     } while( rms_err_too_high );
 *   }
 *   if( should_sync_fields ) synchronize_fields(fields,grid);
 *
 * To clean up the field solver, the following sequence is suggested:
 *   delete_materials(&materials_list);
 *   delete_materials_coefficients(&materials_coefficients);
 *   delete_fields(&fields);
 *****************************************************************************/

typedef struct _field_t {
  float ex, ey, ez;                /* Electric field */
  float cbx, cby, cbz;             /* Magnetic field */
  float tcax, tcay, tcaz;          /* TCA fields */
  float jfx, jfy, jfz;             /* Free current */
  material_id ematx, ematy, ematz; /* Material at edge centers */
  material_id fmatx, fmaty, fmatz; /* Material at face centers */
  material_id nmat, cmat;          /* Material at cell centers and nodes */
  float rhof, rhob;                /* Free and bound charge density */
  float div_e_err, div_b_err;      /* Divergence errors */
} field_t;

typedef struct _hydro_t {
  float rho;           /* Charge density         => < q f > */
  float jx, jy, jz;    /* Current density        => < q v_i f > */
  float ke;            /* Kinetic energy density => < m c^2 (gamma-1) f > */
  float px, py, pz;    /* Momentum density       => < p_i f > */
  float txx, tyy, tzz; /* Stress diagonal        => < p_i v_j f >, i==j */
  float tyz, tzx, txy; /* Stress off-diagonal    => < p_i v_j f >, i!=j */
  float pad0, pad1;    /* 16-byte align the structure */
} hydro_t;

typedef struct _interpolator_t {
  float ex, dexdy, dexdz, d2exdydz;
  float ey, deydz, deydx, d2eydzdx;
  float ez, dezdx, dezdy, d2ezdxdy;
  float cbx, dcbxdx;
  float cby, dcbydy;
  float cbz, dcbzdz;
  float pad0, pad1; /* 16-byte align the structure */
} interpolator_t;

typedef struct _accumulator_t {
  float jx[4]; /* jx0@(0,-1,-1),jx1@(0,1,-1),jx2@(0,-1,1),jx3@(0,1,1) */
  float jy[4]; /* jy0@(-1,0,-1),jy1@(-1,0,1),jy2@(1,0,-1),jy3@(1,0,1) */
  float jz[4]; /* jz0@(-1,-1,0),jz1@(1,-1,0),jz2@(-1,1,0),jz3@(1,1,0) */
} accumulator_t;

BEGIN_C_DECLS

/* In structors.c */

extern field_t * ALIGNED new_field( const grid_t * RESTRICT g );
extern void delete_field( field_t ** ALIGNED f );
extern void clear_field( field_t * RESTRICT ALIGNED f,
		         const grid_t * RESTRICT g );
extern void clear_jf( field_t * RESTRICT ALIGNED f,
		      const grid_t * RESTRICT g );
extern void clear_rhof( field_t * RESTRICT ALIGNED f,
			const grid_t * RESTRICT g );

extern hydro_t * ALIGNED new_hydro( const grid_t * RESTRICT g );
extern void delete_hydro( hydro_t ** ALIGNED h );
extern void clear_hydro( hydro_t * RESTRICT ALIGNED h,
		         const grid_t * RESTRICT g );

extern interpolator_t * ALIGNED new_interpolator( const grid_t * RESTRICT g );
extern void delete_interpolator( interpolator_t ** ALIGNED fi );
extern void clear_interpolator( interpolator_t * RESTRICT ALIGNED fi,
			        const grid_t * RESTRICT g );

extern accumulator_t * ALIGNED new_accumulator( const grid_t * RESTRICT g );
extern void delete_accumulator( accumulator_t ** ALIGNED a );
extern void clear_accumulator( accumulator_t * RESTRICT ALIGNED a,
			       const grid_t * RESTRICT g );

/* In advance_e.c */

extern void compute_curl_b( field_t * RESTRICT ALIGNED f,
			    const material_coefficient_t * RESTRICT ALIGNED m,
			    const grid_t * RESTRICT g );
extern void advance_e( field_t * RESTRICT ALIGNED f,
		       const material_coefficient_t * RESTRICT ALIGNED m,
		       const grid_t * RESTRICT g );

/* In advance_b.c */

extern void advance_b( field_t * RESTRICT ALIGNED f,
		       const grid_t * RESTRICT g,
		       float frac );

/* In energy_f.c */

extern void energy_f( double * RESTRICT energy,
                      const field_t * RESTRICT ALIGNED f,
                      const material_coefficient_t * RESTRICT ALIGNED m,
                      const grid_t * RESTRICT g );

/* In div_e.c */

extern void compute_rhob( field_t * RESTRICT ALIGNED f,
                          const material_coefficient_t * RESTRICT ALIGNED m,
                          const grid_t * RESTRICT g ); 
extern void
compute_div_e_err( field_t * RESTRICT ALIGNED f,
                   const material_coefficient_t * RESTRICT ALIGNED m,
                   const grid_t * RESTRICT g );
extern double compute_rms_div_e_err( field_t * RESTRICT ALIGNED f,
                                     const grid_t * RESTRICT g );
extern void clean_div_e( field_t * RESTRICT ALIGNED f,
			 const material_coefficient_t * RESTRICT ALIGNED m,
			 const grid_t * RESTRICT g );

/* In div_b.c */

extern void compute_div_b_err( field_t * RESTRICT ALIGNED f,
                               const grid_t * RESTRICT g );
extern double compute_rms_div_b_err( field_t * RESTRICT ALIGNED f,
                                     const grid_t * RESTRICT g );
extern void clean_div_b( field_t * RESTRICT ALIGNED f,
			 const grid_t * RESTRICT g );

/* In remote.c */

extern double synchronize_tang_e_norm_b( field_t * RESTRICT ALIGNED f,
					 const grid_t * RESTRICT g );
extern void synchronize_jf( field_t * RESTRICT ALIGNED f,
			    const grid_t * RESTRICT g );
extern void synchronize_rhof( field_t * RESTRICT ALIGNED f,
                              const grid_t * RESTRICT g );
extern void synchronize_rhob( field_t * RESTRICT ALIGNED f,
                              const grid_t * RESTRICT g );
extern void synchronize_hydro( hydro_t * RESTRICT ALIGNED hydro,
			       const grid_t * RESTRICT g );

/* In interpolator.c */

extern void load_interpolator( interpolator_t * RESTRICT ALIGNED fi,
                               const field_t * RESTRICT ALIGNED f,
                               const grid_t * RESTRICT g );

/* In accumulator.c */

extern void unload_accumulator( field_t * RESTRICT ALIGNED f, 
                                const accumulator_t * RESTRICT ALIGNED a,
                                const grid_t * RESTRICT g );

/* field module INTERNAL USE ONLY ********************************************/

/* In local.c */

extern void local_ghost_tang_b( field_t * RESTRICT ALIGNED f,
				const grid_t * RESTRICT g );
extern void local_ghost_norm_e( field_t * RESTRICT ALIGNED f,
                                const grid_t * RESTRICT g );
extern void local_ghost_div_b( field_t * RESTRICT ALIGNED f,
			       const grid_t * RESTRICT g );

extern void local_adjust_tang_e( field_t * RESTRICT ALIGNED f,
                                        const grid_t * RESTRICT g );
extern void local_adjust_div_e( field_t * RESTRICT ALIGNED f,
                                const grid_t * RESTRICT g );
extern void local_adjust_norm_b( field_t * RESTRICT ALIGNED f,
                                 const grid_t * RESTRICT g );

extern void local_adjust_jf( field_t * RESTRICT ALIGNED f,
                             const grid_t * RESTRICT g );
extern void local_adjust_rhof( field_t * RESTRICT ALIGNED f,
                               const grid_t * RESTRICT g );
extern void local_adjust_rhob( field_t * RESTRICT ALIGNED f,
                               const grid_t * RESTRICT g );
extern void local_adjust_hydro( hydro_t * RESTRICT ALIGNED h,
				const grid_t * RESTRICT g );

/* In remote.c */

extern void begin_remote_ghost_tang_b( field_t * RESTRICT ALIGNED f,
				       const grid_t * RESTRICT g );
extern void end_remote_ghost_tang_b( field_t * RESTRICT ALIGNED f,
				     const grid_t * RESTRICT g );
extern void begin_remote_ghost_norm_e( field_t * RESTRICT ALIGNED f,
				       const grid_t * RESTRICT g );
extern void end_remote_ghost_norm_e( field_t * RESTRICT ALIGNED f,
				     const grid_t * RESTRICT g );
extern void begin_remote_ghost_div_b( field_t * RESTRICT ALIGNED f,
				      const grid_t * RESTRICT g );
extern void end_remote_ghost_div_b( field_t * RESTRICT ALIGNED f,
				    const grid_t * RESTRICT g );

END_C_DECLS

#endif
