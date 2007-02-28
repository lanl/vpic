/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version (uses algorithms heavily revised from
 *                    earlier V4PIC versions)
 *
 */

/*****************************************************************************
 * The theory behind the Marder correction is that the Ampere and Faraday
 * equations can be modified as follows:
 *   pB/pt = -curl E    --> pB/pt = -curl E     + alpha grad div B
 *   pD/pt = curl H - J --> pD/pt =  curl H - J + alpha grad ( div D - rho )
 * Taking the divergence of the modified equations yield:
 *   p(div B)/pt     = alpha laplacian div B
 *   p(div D-rho)/pt = alpha laplacian ( div D - rho )
 * Since these are sourceless diffusion equation, asymptotically,
 *   div B       --> 0 
 *   div D - rho --> 0
 * In particular, Fourier transforming div B in space shows that a given mode
 * decays as exp(-alpha k^2 t). The diffusion coefficient alpha controls how
 * fast the divergence goes to zero. Note the long wavelength modes decay more
 * slowly and k=0 does not decay at all. This is not a problem though owing to
 * properties of the divergence operator:
 *   div B @ k=0       -> integral over all space of div B -> 0
 *   div D - rho @ k=0 -> integral over all space of div D - rho
 *                     -> -net charge in universe -> 0
 *
 * Since div B and div D-rho is ideally zero, the modified equations do not
 * change _any_ physics. Further, if for any reason a non-zero div B or
 * (div D - rho) occurs, the above modification will drive the error back to
 * zero.
 *   
 * To understand how use this in a simulation, consider the standard field
 * update equations for Bx on a Yee mesh without the additional term:
 *   cBx(1/2) = cBx(-1/2) - (c*dt)(curl E)_x
 * Because of finite precision arithmetic, cBx(1/2) can be off with a relative
 * error on order the machine's floating point precision, eps (~1.2e-7 for
 * IEEE single precision). Over many time steps, these errors accumulate. The
 * accumulation process can be thought of as a random walk with a RMS step
 * size on the order of ~0.5 eps |cBx|. Thus, over a large number of time
 * steps Nt, the PDF of the error in cBx for an arbitrary grid point will be
 * closely approximated by a Gaussian with zero mean and standard deviation
 * ~0.5 eps |cBx| sqrt(Nt). The same holds true for cBy and cBz.
 * 
 * If it is assumed that the errors between different grid points are
 * uncorrelated (a _very_ accurate assumption except for very specially
 * prepared field configurations), then the power in various spectral modes of
 * div B on a periodic mesh can be shown to be:
 *   |div cB(kx,ky,kz)_unclean|^2 ~
 *     [eps^2 |cB|^2 Nt/(Nx Ny Nz)][ (sin(pi kx/Nx)/dx)^2 +
 *                                   (sin(pi ky/Ny)/dy)^2 +
 *                                   (sin(pi kz/Nz)/dz)^2 ]
 * To reduce this error accumulation, the grad div B term is applied using
 * forward differencing in time (this is the usual Marder pass ... strictly
 * local operations, easy and efficient to implement in parallel):
 *   cBx(1/2)_clean = cBx(1/2)_unclean + 
 *       alpha dt grad div cBx(1/2)_unclean
 * The power in various modes of cBx(1/2)_clean can be shown to be:
 *  |div cB(kx,ky,kz)_clean|^2 ~
 *     |div cB(kx,ky,kz)_unclean|^2 
 *       { 1 - (4*alpha*dt/dg^2) [ (dg sin(pi kx/Nx)/dx)^2 +
 *                                 (dg sin(pi ky/Ny)/dy)^2 +
 *                                 (dg sin(pi kz/Nz)/dz)^2 ] }^2
 * where dg^-2 = dx^-2 + dy^-2 + dz^-2.
 *
 * Inspecting the above, if 0 <= alpha dt < dg^2/2, then no component of
 * div cB(kx,ky,kz) grows and the divergence cleaning pass is numerically
 * stable. Note: This is the same stability criterion as the forward
 * differenced diffusion equation.
 * 
 * If alpha dt = dg^2/4, then shortest wavelength component of div cB will be
 * zeroed. Since this is where most of the divergence errors are located, this
 * is a relatively good choice.
 * 
 * If we want to minimize the total RMS divergence error, it can be shown
 * (using Parseval's theorem) that the best choice of alpha dt on large cubic
 * periodic meshes is:
 *   alpha dt ~ 0.388888889 dg^2
 * This value is pretty close to optimal on other meshes also. Using this
 * value will take the total RMS divergence error to ~0.304 of the original
 * value. 
 * 
 * If we assume future contributions to the divergence error are uncorrelated
 * with previous contributions (a very accurate assumption) and we are only
 * going to clean every Nc time steps, then the maximum relative RMS
 * divergence error, div_max will obey the following relation asymptotically:
 *   div_max^2 = (0.304 div_max)^2    +    0.25 eps^2 Nc
 *                      |                        |
 *                      |                        |
 *       Error left over from previous clean     |
 *                        Error accumulated since previous clean
 * Solving for Nc yields:
 *   Nc ~ 3.63 (div_max/eps)^2
 *
 * Example:
 *   For div_max ~ 1e-6 in single precision, a divergence clean
 *   for cB should be done every 255 time steps.
 *
 * For the clean_div_e, there are two additional considerations. Since
 * advance_e uses exponential differencing in time, the diffusion pass should
 * be modified to be consistent with advance_e. If divergence cleaning were
 * done as part of the time step:
 *   E_new = decay E_old + drive [TCA_new - (dt/eps0) J +
 *                                (dt/eps0)alpha grad(div eps0 epsr E - rho)]
 * Extracting out the divergence cleaning correction yields:
 *   E_clean = E_unclean + drive alpha dt grad (div epsr E - rho/eps0)
 * Second, in conductive medium, the total charge is not known (conduction
 * charge density is not directly computed). Thus, in these regions, the
 * diveregence error is considered to be zero. This gives the complete
 * modified Marder pass:
 *  E_clean = E_unclean +
 *            drive alpha dt grad nonconductive (div epsr E - rho/eps0)
 *****************************************************************************/

#include <math.h> /* For sqrt */
#include <field.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

/*****************************************************************************
 * compute_rhob applies the following difference equation
 *   rho_b = eps0 nonconductive ( div epsr E - rho_f/eps0 )
 * rho_b is not computed correctly on absorbing boundary surfaces but it does
 * not matter as divergence cleaning is not applied there. The structure of
 * this routine is identical to compute_div_e_err, so errors in one will 
 * likely be replicated in the other.
 *****************************************************************************/

void compute_rhob( field_t * RESTRICT ALIGNED f,
                   const material_coefficient_t * RESTRICT ALIGNED m,
                   const grid_t * RESTRICT g ) {
  float px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) { ERROR(("Bad field"));                 return; }
  if( m==NULL ) { ERROR(("Bad material coefficients")); return; }
  if( g==NULL ) { ERROR(("Bad grid"));                  return; }
  
  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? g->eps0/g->dx : 0;
  py = (ny>1) ? g->eps0/g->dy : 0;
  pz = (nz>1) ? g->eps0/g->dz : 0;

# define UPDATE_RHO_B() f0->rhob = m[f0->nmat].nonconductive*      \
    ( px*( m[f0->ematx].epsx*f0->ex - m[fx->ematx].epsx*fx->ex ) + \
      py*( m[f0->ematy].epsy*f0->ey - m[fy->ematy].epsy*fy->ey ) + \
      pz*( m[f0->ematz].epsz*f0->ez - m[fz->ematz].epsz*fz->ez ) - \
      f0->rhof )

  /* Begin setting normal e ghosts */
  begin_remote_ghost_norm_e( f, g );
  local_ghost_norm_e( f, g );

  /* Compute the divergence error in interior */
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(2,y,  z);
      fx = &f(1,y,  z);
      fy = &f(2,y-1,z);
      fz = &f(2,y,  z-1);
      for( x=2; x<=nx; x++ ) {
	UPDATE_RHO_B();
	f0++;
	fx++;
	fy++;
	fz++;
      }
    }
  }
  
  /* Finish setting normal E ghosts */
  end_remote_ghost_norm_e( f, g );

  /* Compute divergence error in exterior */

  /* z faces, x edges, y edges and all corners */
  for( y=1; y<=ny+1; y++ ) {
    f0 = &f(1,y,  1);
    fx = &f(0,y,  1);
    fy = &f(1,y-1,1);
    fz = &f(1,y,  0);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_RHO_B();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }
  for( y=1; y<=ny+1; y++ ) {
    f0 = &f(1,y,  nz+1);
    fx = &f(0,y,  nz+1);
    fy = &f(1,y-1,nz+1);
    fz = &f(1,y,  nz);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_RHO_B();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }
 
  /* y faces, z edges */
  for( z=2; z<=nz; z++ ) {
    f0 = &f(1,1,z);
    fx = &f(0,1,z);
    fy = &f(1,0,z);
    fz = &f(1,1,z-1);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_RHO_B();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }
  for( z=2; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fx = &f(0,ny+1,z);
    fy = &f(1,ny,  z);
    fz = &f(1,ny+1,z-1);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_RHO_B();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }

  /* x faces */
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(1,y,  z);
      fx = &f(0,y,  z);
      fy = &f(1,y-1,z);
      fz = &f(1,y,  z-1);
      UPDATE_RHO_B();
      f0 = &f(nx+1,y,  z);
      fx = &f(nx,  y,  z);
      fy = &f(nx+1,y-1,z);
      fz = &f(nx+1,y,  z-1);
      UPDATE_RHO_B();
    }
  }

  local_adjust_rhob(f,g);
}

/*****************************************************************************
 * compute_div_e_error applies the following difference equation:
 *   div_e_err = nonconductive [ ( div eps_r E ) - ( rho_f + rho_b )/eps_0 ]
 * Note: the structure of this function is identical to compute_rhob so
 * errors in one will likely be replicated in the other.
 *****************************************************************************/

void compute_div_e_err( field_t * RESTRICT ALIGNED f,
                        const material_coefficient_t * RESTRICT ALIGNED m,
                        const grid_t * RESTRICT g ) {
  float px, py, pz, cj;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) { ERROR(("Bad field"));                 return; }
  if( m==NULL ) { ERROR(("Bad material coefficients")); return; }
  if( g==NULL ) { ERROR(("Bad grid"));                  return; }

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? 1./g->dx : 0;
  py = (ny>1) ? 1./g->dy : 0;
  pz = (nz>1) ? 1./g->dz : 0;
  cj = 1./g->eps0;

# define UPDATE_DERR_E() f0->div_e_err = m[f0->nmat].nonconductive * \
    ( px*( m[f0->ematx].epsx*f0->ex - m[fx->ematx].epsx*fx->ex ) +   \
      py*( m[f0->ematy].epsy*f0->ey - m[fy->ematy].epsy*fy->ey ) +   \
      pz*( m[f0->ematz].epsz*f0->ez - m[fz->ematz].epsz*fz->ez ) -   \
      cj*( f0->rhof + f0->rhob ) )

  /* Begin setting normal e ghosts */
  begin_remote_ghost_norm_e( f, g );
  local_ghost_norm_e( f, g );

  /* Compute the divergence error in interior */
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(2,y,  z);
      fx = &f(1,y,  z);
      fy = &f(2,y-1,z);
      fz = &f(2,y,  z-1);
      for( x=2; x<=nx; x++ ) {
	UPDATE_DERR_E();
	f0++;
	fx++;
	fy++;
	fz++;
      }
    }
  }

  /* Finish setting normal E ghosts */
  end_remote_ghost_norm_e( f, g );

  /* Compute divergence error in exterior */

  /* z faces, x edges, y edges and all corners */
  for( y=1; y<=ny+1; y++ ) {
    f0 = &f(1,y,  1);
    fx = &f(0,y,  1);
    fy = &f(1,y-1,1);
    fz = &f(1,y,  0);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_DERR_E();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }
  for( y=1; y<=ny+1; y++ ) {
    f0 = &f(1,y,  nz+1);
    fx = &f(0,y,  nz+1);
    fy = &f(1,y-1,nz+1);
    fz = &f(1,y,  nz);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_DERR_E();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }

  /* y faces, z edges */
  for( z=2; z<=nz; z++ ) {
    f0 = &f(1,1,z);
    fx = &f(0,1,z);
    fy = &f(1,0,z);
    fz = &f(1,1,z-1);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_DERR_E();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }
  for( z=2; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fx = &f(0,ny+1,z);
    fy = &f(1,ny,  z);
    fz = &f(1,ny+1,z-1);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_DERR_E();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }

  /* x faces */
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(1,y,  z);
      fx = &f(0,y,  z);
      fy = &f(1,y-1,z);
      fz = &f(1,y,  z-1);
      UPDATE_DERR_E();
      f0 = &f(nx+1,y,  z);
      fx = &f(nx,  y,  z);
      fy = &f(nx+1,y-1,z);
      fz = &f(nx+1,y,  z-1);
      UPDATE_DERR_E();
    }
  }

  local_adjust_div_e(f,g);
}

/*****************************************************************************
 * compute_rms_div_e_err returns
 *   eps0 sqrt( Integral |div_e_err|^2 / Volume )
 * This has units of (electric charge / volume). The integrals are done over
 * all the domains. The volume is the total volume of all domains. Every
 * processor gets the same value.
 *****************************************************************************/

double compute_rms_div_e_err( field_t * RESTRICT ALIGNED f,
                              const grid_t * RESTRICT g ) {
  double pyz, pz, err, local[2], global[2];
  int x, y, z, nx, ny, nz;
  field_t *f0;

  if( f==NULL ) { ERROR(("Bad field")); return -1; }
  if( g==NULL ) { ERROR(("Bad grid"));  return -1; }
  
  nx = g->nx;
  ny = g->ny;
  nz = g->nz;

  err = 0;
  for( z=1; z<=nz+1; z++ ) {
    pz = (z==1 || z==nz+1) ? 0.5 : 1;
    for( y=1; y<=ny+1; y++ ) {
      pyz = pz*( (y==1 || y==ny+1) ? 0.5 : 1 );
      f0 = &f(1,y,z);
      err += 0.5*pyz*f0->div_e_err*f0->div_e_err;
      f0++;
      for( x=2; x<=nx; x++ ) {
        err += pyz*f0->div_e_err*f0->div_e_err;
        f0++;
      }
      err += 0.5*pyz*f0->div_e_err*f0->div_e_err;
    }
  }

  local[0] = err*g->dx*g->dy*g->dz;
  local[1] = g->nx*g->ny*g->nz*g->dx*g->dy*g->dz;
  mp_allsum_d( local, global, 2, g->mp );
  return g->eps0*sqrt(global[0]/global[1]);
}

/*****************************************************************************
 * clean_div_e applies the following difference equation:
 *   E_new = E_old + drive alpha dt grad div_e_err
 *****************************************************************************/

void clean_div_e( field_t * RESTRICT ALIGNED f,
		  const material_coefficient_t * RESTRICT ALIGNED m,
		  const grid_t * RESTRICT g ) {
  float alphadt, px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) { ERROR(("Bad field"));                 return; }
  if( m==NULL ) { ERROR(("Bad material coefficients")); return; }
  if( g==NULL ) { ERROR(("Bad grid"));                  return; }

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? 1./g->dx : 0;
  py = (ny>1) ? 1./g->dy : 0;
  pz = (nz>1) ? 1./g->dz : 0;
  alphadt = 0.3888889/( px*px + py*py + pz*pz );
  px *= alphadt;
  py *= alphadt;
  pz *= alphadt;
  
# define MARDER_EX() \
    f0->ex += m[f0->ematx].drivex*px*(fx->div_e_err-f0->div_e_err)
# define MARDER_EY() \
    f0->ey += m[f0->ematy].drivey*py*(fy->div_e_err-f0->div_e_err);
# define MARDER_EZ() \
    f0->ez += m[f0->ematz].drivez*pz*(fz->div_e_err-f0->div_e_err);

  /* Do majority in a single pass */
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(1,y,  z);
      fx = &f(2,y,  z);
      fy = &f(1,y+1,z);
      fz = &f(1,y,z+1);
      for( x=1; x<=nx; x++ ) {
	MARDER_EX();
	MARDER_EY();
	MARDER_EZ();
	f0++;
	fx++;
	fy++;
	fz++;
      }
    }
  }

  /* Do left over ex */
  for( y=1; y<=ny+1; y++ ) {
    f0 = &f(1,y,nz+1);
    fx = &f(2,y,nz+1);
    for( x=1; x<=nx; x++ ) {
      MARDER_EX();
      f0++;
      fx++;
    }
  }
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fx = &f(2,ny+1,z);
    for( x=1; x<=nx; x++ ) {
      MARDER_EX();
      f0++;
      fx++;
    }
  }

  /* Do left over ey */
  for( z=1; z<=nz+1; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(nx+1,y,  z);
      fy = &f(nx+1,y+1,z);
      MARDER_EY();
    }
  }
  for( y=1; y<=ny; y++ ) {
    f0 = &f(1,y,  nz+1);
    fy = &f(1,y+1,nz+1);
    for( x=1; x<=nx; x++ ) {
      MARDER_EY();
      f0++;
      fy++;
    }
  }

  /* Do left over ez */
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fz = &f(1,ny+1,z+1);
    for( x=1; x<=nx+1; x++ ) {
      MARDER_EZ();
      f0++;
      fz++;
    }
  }
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(nx+1,y,z);
      fz = &f(nx+1,y,z+1);
      MARDER_EZ();
    }
  }

  local_adjust_tang_e(f,g);
}
