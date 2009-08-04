#define IN_spa
#include "spa_private.h"

// Accumulate particle to rhob with locally adjusted accumulation.
// FIXME: THIS FUNCTION DOESN'T BELONG HERE ANYMORE!

// FIXME: BOUNDARY_P SHOULD BATCH ALL RHOB ACCUMULATIONS TOGETHER.
// WOULD BE LOWER OVERHEAD

#if 0
using namespace v4;

// Note: If part of the body of accumulate_rhob, under the hood
// there is a checked for initialization that occurs everytime
// accumulate_rhob is called!

static const v4float ax[4] = { v4float(1,1,1,1), v4float(2,1,2,1),
                               v4float(1,2,1,2), v4float(2,2,2,2) };
static const v4float ay[4] = { v4float(1,1,1,1), v4float(2,2,1,1),
                               v4float(1,1,2,2), v4float(2,2,2,2) };
#endif

void
accumulate_rhob( field_t          * RESTRICT ALIGNED(128) f,
                 const particle_t * RESTRICT ALIGNED(32)  p,
                 const grid_t     * RESTRICT              g ) {
# if 1

  // See note in rhof for why this variant is used.
  float w0 = p->dx, w1 = p->dy, w2, w3, w4, w5, w6, w7, dz = p->dz;
  int v = p->i, x, y, z, sy = g->sy, sz = g->sz;
  w7 = p->q*g->r8V;

  // Compute the trilinear weights
  // See note in rhof for why FMA and FNMS are done this way.

# define FMA( x,y,z) ((z)+(x)*(y))
# define FNMS(x,y,z) ((z)-(x)*(y))
  w6=FNMS(w0,w7,w7);                    // q(1-dx)
  w7=FMA( w0,w7,w7);                    // q(1+dx)
  w4=FNMS(w1,w6,w6); w5=FNMS(w1,w7,w7); // q(1-dx)(1-dy), q(1+dx)(1-dy)
  w6=FMA( w1,w6,w6); w7=FMA( w1,w7,w7); // q(1-dx)(1+dy), q(1+dx)(1+dy)
  w0=FNMS(dz,w4,w4); w1=FNMS(dz,w5,w5); w2=FNMS(dz,w6,w6); w3=FNMS(dz,w7,w7);
  w4=FMA( dz,w4,w4); w5=FMA( dz,w5,w5); w6=FMA( dz,w6,w6); w7=FMA( dz,w7,w7);
# undef FNMS
# undef FMA

  // Adjust the weights for a corrected local accumulation of rhob.
  // See note in synchronize_rho why we must do this for rhob and not
  // for rhof.

  x  = v;    z = x/sz;
  if( z==1     ) w0 += w0, w1 += w1, w2 += w2, w3 += w3;
  if( z==g->nz ) w4 += w4, w5 += w5, w6 += w6, w7 += w7;
  x -= sz*z; y = x/sy;
  if( y==1     ) w0 += w0, w1 += w1, w4 += w4, w5 += w5;
  if( y==g->ny ) w2 += w2, w3 += w3, w6 += w6, w7 += w7;
  x -= sy*y;
  if( x==1     ) w0 += w0, w2 += w2, w4 += w4, w6 += w6;
  if( x==g->nx ) w1 += w1, w3 += w3, w5 += w5, w7 += w7;

  // Reduce the particle charge to rhob

  f[v      ].rhob += w0; f[v      +1].rhob += w1;
  f[v   +sy].rhob += w2; f[v   +sy+1].rhob += w3;
  f[v+sz   ].rhob += w4; f[v+sz   +1].rhob += w5;
  f[v+sz+sy].rhob += w6; f[v+sz+sy+1].rhob += w7;

# else

  v4float q, wl, wh, rl, rh;
  int v, sy = g->sy, sz = g->sz;
  int i, j;

  // Gather rhob for this voxel

  v = p->i;
  rl = v4float( f[v      ].rhob, f[v      +1].rhob,
                f[v   +sy].rhob, f[v   +sy+1].rhob);
  rh = v4float( f[v+sz   ].rhob, f[v+sz   +1].rhob,
                f[v+sz+sy].rhob, f[v+sz+sy+1].rhob);

  // Compute the trilinear weights

  load_4x1( &p->dx, wl );
  trilinear( wl, wh );

  // Adjust the weights for a corrected local accumulation of rhob.
  // See note in synchronize_rho why we must do this for rhob and not
  // for rhof.  Why yes, this code snippet is branchless and evil.

  i = v;
  j = i/sz; i -= sz*j;       load_4x1( &ax[(j==1    )?3:0], q ); wl *= q;
  /**/                       load_4x1( &ax[(j==g->nz)?3:0], q ); wh *= q;
  j = i/sy; i -= sy*j;
  j = (j==1) + 2*(j==g->ny); load_4x1( &ay[j], q ); wl *= q; wh *= q;
  i = (i==1) + 2*(i==g->nx); load_4x1( &ax[i], q ); wl *= q; wh *= q;

  // Reduce the particle charge to rhof and scatter the result

  q = v4float( g->r8V*p->q );
  store_4x1_tr( fma(q,wl,rl), &f[v      ].rhob, &f[v      +1].rhob,
                              &f[v   +sy].rhob, &f[v   +sy+1].rhob );
  store_4x1_tr( fma(q,wh,rh), &f[v+sz   ].rhob, &f[v+sz   +1].rhob,
                              &f[v+sz+sy].rhob, &f[v+sz+sy+1].rhob );

# endif
}

// FIXME: ARCHITECTURAL FLAW!  CUSTOM BCS AND SHARED FACES CANNOT
// COEXIST ON THE SAME FACE!  THIS MEANS THAT CUSTOM BOUNDARYS MUST
// REINJECT ALL ABSORBED PARTICLES IN THE SAME DOMAIN!

using namespace v4;

void
boundary_p( species_t     * RESTRICT sp_list,
            field_t       * RESTRICT ALIGNED(128) f,
            accumulator_t * RESTRICT ALIGNED(128) a0,
            const grid_t  * RESTRICT g,
            mt_rng_t      * RESTRICT rng ) {

  // Gives the local mp buffer associated with a local face
  static const int f2b[6]  = { BOUNDARY(-1, 0, 0),
                               BOUNDARY( 0,-1, 0),
                               BOUNDARY( 0, 0,-1),
                               BOUNDARY( 1, 0, 0),
                               BOUNDARY( 0, 1, 0),
                               BOUNDARY( 0, 0, 1) };

  // Gives the remote mp buffer associated with a local face
  static const int f2rb[6] = { BOUNDARY( 1, 0, 0),
                               BOUNDARY( 0, 1, 0),
                               BOUNDARY( 0, 0, 1),
                               BOUNDARY(-1, 0, 0),
                               BOUNDARY( 0,-1, 0),
                               BOUNDARY( 0, 0,-1) };

  // Temporary store for local particle injectors
  // FIXME: Ugly static usage
  static particle_injector_t * RESTRICT ALIGNED(16) ci = NULL;
  static int max_ci = 0;

  mp_handle mp    = g->mp;
  const int rank  = mp_rank(mp);
  const int nproc = mp_nproc(mp);
  /*const*/ int bc[6], shared[6];
  int face, n_send[6], n_recv[6], n_ci;

  // Setup and begin receiving the particle counts

  for( face=0; face<6; face++ ) {
    bc[face] = g->bc[f2b[face]];
    shared[face] = (bc[face]>=0) && (bc[face]<nproc) && (bc[face]!=rank);
    if( shared[face] ) {
      mp_size_recv_buffer( f2b[face], sizeof(int), mp );
      mp_begin_recv( f2b[face], sizeof(int), bc[face], f2rb[face], mp );
    }
  }

  // Load the particle send and local injection buffers
  
  do {
    static const int axis[6]  = { 0, 1, 2,  0,  1,  2 };
    static const float dir[6] = { 1, 1, 1, -1, -1, -1 };

    const int64_t    * RESTRICT ALIGNED(128) neighbor = g->neighbor;
    /**/  boundary_t * RESTRICT              boundary = g->boundary;
    const int                                nb       = g->nb; 
    const int64_t                            rangel   = g->rangel;
    const int64_t                            rangeh   = g->rangeh;
    const int64_t                            rangem   = g->range[nproc];
    /**/  int64_t range[6];

    particle_injector_t * RESTRICT ALIGNED(16) pi_send[6];
    particle_injector_t * RESTRICT ALIGNED(16) pi;
    species_t * RESTRICT sp;
    int i, voxel, np, nm;
    int64_t nn;

    // Presize the send and injection buffers
    //
    // Each buffer is large enough to hold one injector corresponding
    // to every mover in use (worst case, but plausible scenario in
    // beam simulations, is one buffer gets all the movers).
    //
    // Note: We could be several times more efficient in our particle
    // injector buffer sizing here.  Namely, we could create on local
    // injector buffer of nm is size.  All injection for all
    // boundaries would be done here.  The local buffer would then be
    // counted to determine the size of each send buffer.  The local
    // buffer would then move all injectors into the approate send
    // buffers (leaving only the local injectors).  This would require
    // some extra data motion though.  (But would give a more robust
    // implementation against variations in MP implementation.)
    //
    // FIXME: THIS PRESIZING ASSUMES THAT CUSTOM BOUNDARY CONDITIONS
    // INJECT AT MOST ONE PARTICLE PER INCIDENT PARTICLE.  THIS IS
    // USUALLY TRUE (SOME SECONDARY MODELS MAY NOT SATISFY) BUT THIS
    // ARCHITECTURAL FLAW SHOULD BE FIXED.
    
    nm = 0; LIST_FOR_EACH( sp, sp_list ) nm += sp->nm;
    for( face=0; face<6; face++ )
      if( shared[face] ) {
        mp_size_send_buffer(f2b[face], 16+nm*sizeof(particle_injector_t), mp);
        pi_send[face] =
          (particle_injector_t *)(((char *)mp_send_buffer(f2b[face], mp))+16);
        n_send[face]  = 0;
        range[face]   = g->range[bc[face]];
      }
    if( max_ci<nm ) {
      particle_injector_t * new_pi = ci;
      FREE_ALIGNED( new_pi );
      MALLOC_ALIGNED( new_pi, nm, 16 );
      ci     = new_pi;
      max_ci = nm;
    }
    n_ci = 0;

    // For each species, load the movers

    LIST_FOR_EACH( sp, sp_list ) {
      particle_t       * RESTRICT ALIGNED(128) p0 = sp->p;
      particle_mover_t * RESTRICT ALIGNED(16)  pm = sp->pm + sp->nm - 1;
      const int32_t sp_id = sp->id;
      np = sp->np;
      nm = sp->nm;
      
      // Note that particle movers for each species are processed in
      // reverse order.  This allows us to backfill holes in the
      // particle list created by boundary conditions and/or
      // communication.  This assumes particle on the mover list are
      // monotonically increasing.  That is: pm[n].i > pm[n-1].i for
      // n=1...nm-1.  advance_p and inject_particle create movers with
      // property if all aged particle injection occurs after
      // advance_p and before this

      for( ; nm; pm--, nm-- ) {
        i = pm->i;
        voxel = p0[i].i;
        face = voxel & 7;
        voxel >>= 3;
        p0[i].i = voxel;
        nn = neighbor[ 6*voxel + face ];
        
#       define BACKFILL()                       \
        np--;                                   \
        copy_4x1( &p0[i].dx, &p0[np].dx );      \
        copy_4x1( &p0[i].ux, &p0[np].ux )

        // Absorb

        if( nn==absorb_particles ) {
          accumulate_rhob( f, p0+i, g );
          BACKFILL();
          continue;
        }

        // Send to a neighboring node

        if( ((nn>=0) & (nn< rangel)) | ((nn>rangeh) & (nn<=rangem)) ) {
          pi = &pi_send[face][n_send[face]++];
          copy_4x1( &pi->dx,    &p0[i].dx  );
          copy_4x1( &pi->ux,    &p0[i].ux  );
          copy_4x1( &pi->dispx, &pm->dispx );
          (&pi->dx)[axis[face]] = dir[face];
          pi->i                 = nn - range[face];
          pi->sp_id             = sp_id;
          BACKFILL();
          continue;
        }

        // User-defined handling

        // FIXME: Currently, after a particle interacts with a
        // boundary it is removed from the local particle list.  Thus,
        // if a boundary handler does not want a particle destroyed,
        // it is the boundary handler's job to append the destroyed
        // particle to the list of particles to inject.
        //
        // Note that these destructing and creation processes do _not_
        // adjust rhob by default.  Thus, a boundary handler is
        // responsible for insuring that the rhob is updated
        // appropriate for the incident particle it destroys and for
        // any particles it injects as a result too.
        //
        // Since most boundary handlers do local reinjection and are
        // charge neutral, this means most boundary handlers do
        // nothing to rhob.
        //
        // In the future, the boundary handlers should be adjusted to
        // determine whether or not the interacting particle should be
        // removed in addition to keeping the charge conservation
        // books.

        nn = -nn - 3; // Assumes reflective/absorbing are -1, -2
        if( (nn>=0) & (nn<nb) ) {
          n_ci += boundary[nn].handler( boundary[nn].params, p0+i, pm, f, a0,
                                        g, sp, ci+n_ci, rng, face );
          BACKFILL();
          continue;
        }

        // Uh-oh: We fell through

        WARNING(( "Unknown boundary interaction ... using absorption "
                  "(species=%s, rank=%i)", sp->name, rank ));
        accumulate_rhob( f, p0+i, g );
        BACKFILL();

#       undef BACKFILL
      }
      
      sp->np = np;
      sp->nm = 0;
    }

  } while(0);

  // Finish exchanging particle counts and start exchanging actual
  // particles.
  
  // Note: This is wasteful of communications.  A better protocol
  // would fuse the exchange of the counts with the exchange of the
  // messages.  in a slightly more complex protocol.  Hwoever, the MP
  // API prohibits such a model.  Unfortuantely, refining MP is not
  // much help here.  Under the hood on Roadrunner, the DaCS API also
  // prohibits such (specifically, in both, you can't do the
  // equilvanet of a MPI_Getcount to determine how much data you
  // actually received.
  
  for( face=0; face<6; face++ )
    if( shared[face] ) {
      *((int *)mp_send_buffer( f2b[face], mp )) = n_send[face];
      mp_begin_send( f2b[face], sizeof(int), bc[face], f2b[face], mp );
    }

  for( face=0; face<6; face++ )
    if( shared[face] )  {
      mp_end_recv( f2b[face], mp );
      n_recv[face] = *((int *)mp_recv_buffer( f2b[face], mp ));
      mp_size_recv_buffer( f2b[face],
                           16+n_recv[face]*sizeof(particle_injector_t), mp );
      mp_begin_recv( f2b[face], 16+n_recv[face]*sizeof(particle_injector_t),
                     bc[face], f2rb[face], mp );
    }

  for( face=0; face<6; face++ )
    if( shared[face] ) {
      mp_end_send( f2b[face], mp );
      // FIXME: ASSUMES MP WON'T MUCK WITH REST OF SEND BUFFER IF WE
      // DID MORE EFFICIENT MOVER ALLOCATION ABOVE, THIS WOULD BE
      // ROBUSTED AGAINST VAGARIES OF MP IMPLEMENTATIONS
      mp_begin_send( f2b[face], 16+n_send[face]*sizeof(particle_injector_t),
                     bc[face], f2b[face], mp );
    }
  
  // Inject particles as we finish receiving them

  do {
    const particle_injector_t * RESTRICT ALIGNED(16) pi;
    particle_t       * RESTRICT ALIGNED(32) sp_p[ 64]; int sp_np[64];
    particle_mover_t * RESTRICT ALIGNED(32) sp_pm[64]; int sp_nm[64];
    species_t * sp;
    int n;
    
#   ifdef DISABLE_DYNAMIC_RESIZING
    int sp_max_np[64], n_dropped_particles[64];
    int sp_max_nm[64], n_dropped_movers[64];
#   else
    // Resize each species's particle and mover storage to be large
    // enough to guarantee successful injection.  (If we broke down
    // the n_recv[face] by species before sending it, we could be
    // tighter on memory footprint here.)
    
    int max_inj = n_ci;
    for( face=0; face<6; face++ )
      if( shared[face] ) max_inj += n_recv[face];

    LIST_FOR_EACH( sp, sp_list ) {
      particle_mover_t * new_pm;
      particle_t * new_p;
      
      n = sp->np + max_inj;
      if( n>sp->max_np ) {
        n = n + (n>>2) + (n>>4); // Increase by 31.25% (~<"silver
        /**/                     // ratio") to minimize resizes (max
        /**/                     // rate that avoids excessive heap
        /**/                     // fragmentation)
        WARNING(( "Resizing local %s particle storage from %i to %i",
                  sp->name, sp->max_np, n ));
        MALLOC_ALIGNED( new_p, n, 128 );
        COPY( new_p, sp->p, sp->np );
        FREE_ALIGNED( sp->p );
        sp->p = new_p, sp->max_np = n;
      }
      
      n = sp->nm + max_inj;
      if( n>sp->max_nm ) {
        n = n + (n>>2) + (n>>4); // See note above
        WARNING(( "Resizing local %s mover storage from %i to %i",
                  sp->name, sp->max_nm, n ));
        MALLOC_ALIGNED( new_pm, n, 128 );
        COPY( new_pm, sp->pm, sp->nm );
        FREE_ALIGNED( sp->pm );
        sp->pm = new_pm, sp->max_nm = n;
      }
    }
#   endif

    // Inject particles.  We do custom local injection first to
    // increase message overlap opportunities.

    LIST_FOR_EACH( sp, sp_list ) {
      if( sp->id<0 || sp->id>=64 ) ERROR(( "Invalid sp->id" ));
      sp_p[ sp->id]=sp->p,  sp_pm[sp->id]=sp->pm;
      sp_np[sp->id]=sp->np, sp_nm[sp->id]=sp->nm;
#     ifdef DISABLE_DYNAMIC_RESIZING
      sp_max_np[sp->id]=sp->max_np, sp_max_nm[sp->id]=sp->max_nm;
      n_dropped_particles[sp->id]=0, n_dropped_movers[sp->id]=0;
#     endif
    }

    face = 5;
    do {
      particle_t       * RESTRICT ALIGNED(32) p;
      particle_mover_t * RESTRICT ALIGNED(16) pm;
      int np, nm, id;

      face++; if( face==7 ) face = 0;
      if( face==6 ) pi = ci, n = n_ci;
      else if( shared[face] ) {
        mp_end_recv( f2b[face], mp );
        pi = (const particle_injector_t *)
          (((char *)mp_recv_buffer(f2b[face],mp))+16);
        n  = n_recv[face];
      } else continue;
      
      // Reverse order injection is done to reduce thrashing of the
      // particle list (particles are removed reverse order so the
      // overall impact of removal + injection is to keep injected
      // particles in order).
      //
      // WARNING: THIS TRUSTS THAT THE INJECTORS (INCLUDING THOSE
      // RECEIVED FROM OTHER NODES) HAVE VALID PARTICLE IDS.

      pi += n-1;
      for( ; n; pi--, n-- ) {
        id = pi->sp_id;
        p  = sp_p[id];  np = sp_np[id];
        pm = sp_pm[id]; nm = sp_nm[id];
#       ifdef DISABLE_DYNAMIC_RESIZING
        if( np>=sp_max_np[id] ) { n_dropped_particles[id]++; continue; }
#       endif
        copy_4x1(  &p[np].dx,    &pi->dx    );
        copy_4x1(  &p[np].ux,    &pi->ux    );
        sp_np[id] = np+1;
#       ifdef DISABLE_DYNAMIC_RESIZING
        if( nm>=sp_max_nm[id] ) { n_dropped_movers[id]++;    continue; }
#       endif
        copy_4x1( &pm[nm].dispx, &pi->dispx ); pm[nm].i = np;
        sp_nm[id] = nm + move_p( p, pm+nm, a0, g );
      }
    } while(face!=5);

    LIST_FOR_EACH( sp, sp_list ) {
#     ifdef DISABLE_DYNAMIC_RESIZING
      if( n_dropped_particles[sp->id] )
        WARNING(( "Dropped %i particles from species \"%s\".  Use a larger "
                  "local particle allocation in your simulation setup for "
                  "this species on this node.",
                  n_dropped_particles[sp->id], sp->name ));

      if( n_dropped_movers[sp->id] )
        WARNING(( "%i particles were not completed moved to their final "
                  "location this timestep for species \"%s\".  Use a larger "
                  "local particle mover buffer in your simulation setup "
                  "for this species on this node.",
                  n_dropped_movers[sp->id], sp->name ));
#     endif
      sp->np=sp_np[sp->id], sp->nm=sp_nm[sp->id];
    }

  } while(0);
  
  for( face=0; face<6; face++ )
    if( shared[face] )
      mp_end_send( f2b[face], mp );
}
