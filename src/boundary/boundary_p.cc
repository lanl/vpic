#define IN_boundary
#include "boundary_private.h"

// If this is defined particle and mover buffers will not resize dynamically
// (This is the common case for the users)
//#define DISABLE_DYNAMIC_RESIZING

// FIXME: ARCHITECTURAL FLAW!  CUSTOM BCS AND SHARED FACES CANNOT
// COEXIST ON THE SAME FACE!  THIS MEANS THAT CUSTOM BOUNDARYS MUST
// REINJECT ALL ABSORBED PARTICLES IN THE SAME DOMAIN!


// Updated by Scott V. Luedtke, XCP-6, December 6, 2018.
// The mover array is now resized along with the particle array.  The mover
// array is filled during advance_p and is most likely to overflow there, not
// here.  Both arrays will now resize down as well.
// 12/20/18: The mover array is no longer resized with the particle array, as
// this actually uses more RAM than having static mover arrays.  The mover will
// still size up if there are too many incoming particles, but I have not
// encountered this.  Some hard-to-understand bit shifts have been replaced with
// cleaner code that the compiler should have no trouble optimizing.
// Spits out lots of warnings. TODO: Remove warnings after testing.

#ifdef V4_ACCELERATION
using namespace v4;
#endif

#ifndef MIN_NP
#define MIN_NP 128 // Default to 4kb (~1 page worth of memory)
//#define MIN_NP 32768 // 32768 particles is 1 MiB of memory.
#endif


enum { MAX_PBC = 32, MAX_SP = 32 };

void
boundary_p( particle_bc_t       * RESTRICT pbc_list,
            species_t           * RESTRICT sp_list,
            field_array_t       * RESTRICT fa,
            accumulator_array_t * RESTRICT aa ) {

  // Gives the local mp port associated with a local face
  static const int f2b[6]  = { BOUNDARY(-1, 0, 0),
                               BOUNDARY( 0,-1, 0),
                               BOUNDARY( 0, 0,-1),
                               BOUNDARY( 1, 0, 0),
                               BOUNDARY( 0, 1, 0),
                               BOUNDARY( 0, 0, 1) };

  // Gives the remote mp port associated with a local face
  static const int f2rb[6] = { BOUNDARY( 1, 0, 0),
                               BOUNDARY( 0, 1, 0),
                               BOUNDARY( 0, 0, 1),
                               BOUNDARY(-1, 0, 0),
                               BOUNDARY( 0,-1, 0),
                               BOUNDARY( 0, 0,-1) };

  // Gives the axis associated with a local face
  static const int axis[6]  = { 0, 1, 2,  0,  1,  2 };

  // Gives the location of sending face on the receiver
  static const float dir[6] = { 1, 1, 1, -1, -1, -1 };

  // Temporary store for local particle injectors
  // FIXME: Ugly static usage
  static particle_injector_t * RESTRICT ALIGNED(16) ci = NULL;
  static int max_ci = 0;

  int n_send[6], n_recv[6], n_ci;

  species_t * sp;
  int face;

  // Check input args

  if( !sp_list ) return; // Nothing to do if no species
  if( !fa || !aa || sp_list->g!=aa->g || fa->g!=aa->g )
    ERROR(( "Bad args" ));

  // Unpack the particle boundary conditions

  particle_bc_func_t pbc_interact[MAX_PBC];
  void * pbc_params[MAX_PBC];
  const int nb = num_particle_bc( pbc_list );
  if( nb>MAX_PBC ) ERROR(( "Update this to support more particle boundary conditions" ));
  for( particle_bc_t * pbc=pbc_list; pbc; pbc=pbc->next ) {
    pbc_interact[-pbc->id-3] = pbc->interact;
    pbc_params[  -pbc->id-3] = pbc->params;
   }

  // Unpack fields

  field_t * RESTRICT ALIGNED(128) f = fa->f;
  grid_t  * RESTRICT              g = fa->g;

  // Unpack accumulator

  accumulator_t * RESTRICT ALIGNED(128) a0 = aa->a;

  // Unpack the grid

  const int64_t * RESTRICT ALIGNED(128) neighbor = g->neighbor;
  /**/  mp_t    * RESTRICT              mp       = g->mp;
  const int64_t rangel = g->rangel;
  const int64_t rangeh = g->rangeh;
  const int64_t rangem = g->range[world_size];
  /*const*/ int bc[6], shared[6];
  /*const*/ int64_t range[6];
  for( face=0; face<6; face++ ) {
    bc[face] = g->bc[f2b[face]];
    shared[face] = (bc[face]>=0) && (bc[face]<world_size) &&
                   (bc[face]!=world_rank);
    if( shared[face] ) range[face] = g->range[bc[face]];
  }

  // Begin receiving the particle counts

  for( face=0; face<6; face++ )
    if( shared[face] ) {
      mp_size_recv_buffer( mp, f2b[face], sizeof(int) );
      mp_begin_recv( mp, f2b[face], sizeof(int), bc[face], f2rb[face] );
    }

  // Load the particle send and local injection buffers

  do {

    particle_injector_t * RESTRICT ALIGNED(16) pi_send[6];

    // Presize the send and injection buffers
    //
    // Each buffer is large enough to hold one injector corresponding
    // to every mover in use (worst case, but plausible scenario in
    // beam simulations, is one buffer gets all the movers).
    //
    // FIXME: We could be several times more efficient in our particle
    // injector buffer sizing here.  Namely, we could create on local
    // injector buffer of nm is size.  All injection for all
    // boundaries would be done here.  The local buffer would then be
    // counted to determine the size of each send buffer.  The local
    // buffer would then move all injectors into the approate send
    // buffers (leaving only the local injectors).  This would require
    // some extra data motion though.  (But would give a more robust
    // implementation against variations in MP implementation.)
    //
    // FIXME: This presizing assumes that custom boundary conditions
    // inject at most one particle per incident particle.  Currently,
    // the invocation of pbc_interact[*] insures that assumption will
    // be satisfied (if the handlers conform that it).  We should be
    // more flexible though in the future (especially given above the
    // above overalloc).

    int nm = 0; LIST_FOR_EACH( sp, sp_list ) nm += sp->nm;

    for( face=0; face<6; face++ )
      if( shared[face] ) {
        mp_size_send_buffer( mp, f2b[face], 16+nm*sizeof(particle_injector_t) );
        pi_send[face] = (particle_injector_t *)(((char *)mp_send_buffer(mp,f2b[face]))+16);
        n_send[face] = 0;
      }

    if( max_ci<nm ) {
      particle_injector_t * new_ci = ci;
      FREE_ALIGNED( new_ci );
      MALLOC_ALIGNED( new_ci, nm, 16 );
      ci     = new_ci;
      max_ci = nm;
    }
    n_ci = 0;

    // For each species, load the movers

    LIST_FOR_EACH( sp, sp_list ) {
      const float   sp_q  = sp->q;
      const int32_t sp_id = sp->id;

      particle_t * RESTRICT ALIGNED(128) p0 = sp->p;
      int np = sp->np;

      particle_mover_t * RESTRICT ALIGNED(16)  pm = sp->pm + sp->nm - 1;
      nm = sp->nm;

      particle_injector_t * RESTRICT ALIGNED(16) pi;
      int i, voxel;
      int64_t nn;

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

        // Absorb

        if( nn==absorb_particles ) {
          // Ideally, we would batch all rhob accumulations together
          // for efficiency
          accumulate_rhob( f, p0+i, g, sp_q );
          goto backfill;
        }

        // Send to a neighboring node

        if( ((nn>=0) & (nn< rangel)) | ((nn>rangeh) & (nn<=rangem)) ) {
          pi = &pi_send[face][n_send[face]++];
#         ifdef V4_ACCELERATION
          copy_4x1( &pi->dx,    &p0[i].dx  );
          copy_4x1( &pi->ux,    &p0[i].ux  );
          copy_4x1( &pi->dispx, &pm->dispx );
#         else
          pi->dx=p0[i].dx; pi->dy=p0[i].dy; pi->dz=p0[i].dz;
          pi->ux=p0[i].ux; pi->uy=p0[i].uy; pi->uz=p0[i].uz; pi->w=p0[i].w;
          pi->dispx = pm->dispx; pi->dispy = pm->dispy; pi->dispz = pm->dispz;
#         endif
          (&pi->dx)[axis[face]] = dir[face];
          pi->i                 = nn - range[face];
          pi->sp_id             = sp_id;
          goto backfill;
        }

        // User-defined handling

        // After a particle interacts with a boundary it is removed
        // from the local particle list.  Thus, if a boundary handler
        // does not want a particle destroyed,  it is the boundary
        // handler's job to append the destroyed particle to the list
        // of particles to inject.
        //
        // Note that these destruction and creation processes do _not_
        // adjust rhob by default.  Thus, a boundary handler is
        // responsible for insuring that the rhob is updated
        // appropriate for the incident particle it destroys and for
        // any particles it injects as a result too.
        //
        // Since most boundary handlers do local reinjection and are
        // charge neutral, this means most boundary handlers do
        // nothing to rhob.

        nn = -nn - 3; // Assumes reflective/absorbing are -1, -2
        if( (nn>=0) & (nn<nb) ) {
          n_ci += pbc_interact[nn]( pbc_params[nn], sp, p0+i, pm,
                                    ci+n_ci, 1, face );
          goto backfill;
        }

        // Uh-oh: We fell through

        WARNING(( "Unknown boundary interaction ... dropping particle "
                  "(species=%s)", sp->name ));

      backfill:

        np--;
#       ifdef V4_ACCELERATION
        copy_4x1( &p0[i].dx, &p0[np].dx );
        copy_4x1( &p0[i].ux, &p0[np].ux );
#       else
        p0[i] = p0[np];
#       endif

      }

      sp->np = np;
      sp->nm = 0;
    }

  } while(0);

  // Finish exchanging particle counts and start exchanging actual
  // particles.

  // Note: This is wasteful of communications.  A better protocol
  // would fuse the exchange of the counts with the exchange of the
  // messages.  in a slightly more complex protocol.  However, the MP
  // API prohibits such a model.  Unfortuantely, refining MP is not
  // much help here.  Under the hood on Roadrunner, the DaCS API also
  // prohibits such (specifically, in both, you can't do the
  // equilvanet of a MPI_Getcount to determine how much data you
  // actually received.

  for( face=0; face<6; face++ )
    if( shared[face] ) {
      *((int *)mp_send_buffer( mp, f2b[face] )) = n_send[face];
      mp_begin_send( mp, f2b[face], sizeof(int), bc[face], f2b[face] );
    }

  for( face=0; face<6; face++ )
    if( shared[face] )  {
      mp_end_recv( mp, f2b[face] );
      n_recv[face] = *((int *)mp_recv_buffer( mp, f2b[face] ));
      mp_size_recv_buffer( mp, f2b[face],
                           16+n_recv[face]*sizeof(particle_injector_t) );
      mp_begin_recv( mp, f2b[face], 16+n_recv[face]*sizeof(particle_injector_t),
                     bc[face], f2rb[face] );
    }

  for( face=0; face<6; face++ )
    if( shared[face] ) {
      mp_end_send( mp, f2b[face] );
      // FIXME: ASSUMES MP WON'T MUCK WITH REST OF SEND BUFFER. IF WE
      // DID MORE EFFICIENT MOVER ALLOCATION ABOVE, THIS WOULD BE
      // ROBUSTED AGAINST MP IMPLEMENTATION VAGARIES
      mp_begin_send( mp, f2b[face], 16+n_send[face]*sizeof(particle_injector_t),
                     bc[face], f2b[face] );
    }

# ifndef DISABLE_DYNAMIC_RESIZING
  // Resize particle storage to accomodate worst case inject

  do {
    int n, nm;

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
        n += 0.3125*n; // Increase by 31.25% (~<"silver
        /**/                     // ratio") to minimize resizes (max
        /**/                     // rate that avoids excessive heap
        /**/                     // fragmentation)
        //float resize_ratio = (float)n/sp->max_np;
        WARNING(( "Resizing local %s particle storage from %i to %i",
                  sp->name, sp->max_np, n ));
        MALLOC_ALIGNED( new_p, n, 128 );
        COPY( new_p, sp->p, sp->np );
        FREE_ALIGNED( sp->p );
        sp->p = new_p, sp->max_np = n;

        /*nm = sp->max_nm * resize_ratio;
        WARNING(( "Resizing local %s mover storage from %i to %i",
                  sp->name, sp->max_nm, nm ));
        MALLOC_ALIGNED( new_pm, nm, 128 );
        COPY( new_pm, sp->pm, sp->nm );
        FREE_ALIGNED( sp->pm );
        sp->pm = new_pm;
        sp->max_nm = nm;*/
      }
      else if(sp->max_np > MIN_NP && n < sp->max_np>>1)
      {
        n += 0.125*n; // Overallocate by less since this rank is decreasing
        if (n<MIN_NP) n = MIN_NP;
        //float resize_ratio = (float)n/sp->max_np;
        WARNING(( "Resizing (shrinking) local %s particle storage from "
                    "%i to %i", sp->name, sp->max_np, n));
        MALLOC_ALIGNED( new_p, n, 128 );
        COPY( new_p, sp->p, sp->np );
        FREE_ALIGNED( sp->p );
        sp->p = new_p, sp->max_np = n;

        /*nm = sp->max_nm * resize_ratio;
        WARNING(( "Resizing (shrinking) local %s mover storage from "
                    "%i to %i", sp->name, sp->max_nm, nm));
        MALLOC_ALIGNED( new_pm, nm, 128 );
        COPY( new_pm, sp->pm, sp->nm );
        FREE_ALIGNED( sp->pm );
        sp->pm = new_pm, sp->max_nm = nm;*/
      }

      // Feasibly, a vacuum-filled rank may receive a shock and need more movers
      // than available from MIN_NP
      nm = sp->nm + max_inj;
      if( nm>sp->max_nm ) {
        nm += 0.3125*nm; // See note above
        //float resize_ratio = (float)nm/sp->max_nm;
        WARNING(( "This happened.  Resizing local %s mover storage from "
                    "%i to %i based on not enough movers",
                  sp->name, sp->max_nm, nm ));
        MALLOC_ALIGNED( new_pm, nm, 128 );
        COPY( new_pm, sp->pm, sp->nm );
        FREE_ALIGNED( sp->pm );
        sp->pm = new_pm;
        sp->max_nm = nm;

        /*n = sp->max_np * resize_ratio;
        WARNING(( "Resizing local %s particle storage from %i to %i",
                  sp->name, sp->max_np, n ));
        MALLOC_ALIGNED( new_p, n, 128 );
        COPY( new_p, sp->p, sp->np );
        FREE_ALIGNED( sp->p );
        sp->p = new_p, sp->max_np = n;*/
      }
    }
  } while(0);
# endif

  do {

    // Unpack the species list for random acesss

    particle_t       * RESTRICT ALIGNED(32) sp_p[ MAX_SP];
    particle_mover_t * RESTRICT ALIGNED(32) sp_pm[MAX_SP];
    float sp_q[MAX_SP];
    int sp_np[MAX_SP];
    int sp_nm[MAX_SP];

#   ifdef DISABLE_DYNAMIC_RESIZING
    int sp_max_np[64], n_dropped_particles[64];
    int sp_max_nm[64], n_dropped_movers[64];
#   endif

    if( num_species( sp_list ) > MAX_SP )
      ERROR(( "Update this to support more species" ));
    LIST_FOR_EACH( sp, sp_list ) {
      sp_p[  sp->id ] = sp->p;
      sp_pm[ sp->id ] = sp->pm;
      sp_q[  sp->id ] = sp->q;
      sp_np[ sp->id ] = sp->np;
      sp_nm[ sp->id ] = sp->nm;
#     ifdef DISABLE_DYNAMIC_RESIZING
      sp_max_np[sp->id]=sp->max_np; n_dropped_particles[sp->id]=0;
      sp_max_nm[sp->id]=sp->max_nm; n_dropped_movers[sp->id]=0;
#     endif
    }

    // Inject particles.  We do custom local injection first to
    // increase message overlap opportunities.

    face = 5;
    do {
      /**/  particle_t          * RESTRICT ALIGNED(32) p;
      /**/  particle_mover_t    * RESTRICT ALIGNED(16) pm;
      const particle_injector_t * RESTRICT ALIGNED(16) pi;
      int np, nm, n, id;

      face++; if( face==7 ) face = 0;
      if( face==6 ) pi = ci, n = n_ci;
      else if( shared[face] ) {
        mp_end_recv( mp, f2b[face] );
        pi = (const particle_injector_t *)
          (((char *)mp_recv_buffer(mp,f2b[face]))+16);
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
#       ifdef V4_ACCELERATION
        copy_4x1(  &p[np].dx,    &pi->dx    );
        copy_4x1(  &p[np].ux,    &pi->ux    );
#       else
        p[np].dx=pi->dx; p[np].dy=pi->dy; p[np].dz=pi->dz; p[np].i=pi->i;
        p[np].ux=pi->ux; p[np].uy=pi->uy; p[np].uz=pi->uz; p[np].w=pi->w;
#       endif
        sp_np[id] = np+1;

#       ifdef DISABLE_DYNAMIC_RESIZING
        if( nm>=sp_max_nm[id] ) { n_dropped_movers[id]++;    continue; }
#       endif
#       ifdef V4_ACCELERATION
        copy_4x1( &pm[nm].dispx, &pi->dispx );
        pm[nm].i = np;
#       else
        pm[nm].dispx=pi->dispx; pm[nm].dispy=pi->dispy; pm[nm].dispz=pi->dispz;
        pm[nm].i=np;
#       endif
        sp_nm[id] = nm + move_p( p, pm+nm, a0, g, sp_q[id] );
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
      sp->np=sp_np[sp->id];
      sp->nm=sp_nm[sp->id];
    }

  } while(0);

  for( face=0; face<6; face++ )
    if( shared[face] ) mp_end_send(mp,f2b[face]);
}
