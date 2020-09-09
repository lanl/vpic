#define IN_boundary

#include <iostream> // TODO: delete
#include "boundary_private.h"

// If this is defined particle and mover buffers will not resize dynamically.
// This is the common case for the users.

// #define DISABLE_DYNAMIC_RESIZING

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

#ifdef V8_ACCELERATION
using namespace v8;
#endif

#ifndef MIN_NP
#define MIN_NP 128 // Default to 4kb (~1 page worth of memory)
//#define MIN_NP 32768 // 32768 particles is 1 MiB of memory.
#endif

enum { MAX_PBC = 32, MAX_SP = 32 };

// This is the AoS implementation.

void
boundary_p( particle_bc_t       * RESTRICT pbc_list,
            species_t           * RESTRICT sp_list,
            field_array_t       * RESTRICT fa,
            accumulator_array_t * RESTRICT aa )
{
  // Gives the local mp port associated with a local face.
  static const int f2b[6]  = { BOUNDARY(-1, 0, 0),
                               BOUNDARY( 0,-1, 0),
                               BOUNDARY( 0, 0,-1),
                               BOUNDARY( 1, 0, 0),
                               BOUNDARY( 0, 1, 0),
                               BOUNDARY( 0, 0, 1) };

  // Gives the remote mp port associated with a local face.
  static const int f2rb[6] = { BOUNDARY( 1, 0, 0),
                               BOUNDARY( 0, 1, 0),
                               BOUNDARY( 0, 0, 1),
                               BOUNDARY(-1, 0, 0),
                               BOUNDARY( 0,-1, 0),
                               BOUNDARY( 0, 0,-1) };

  // Gives the axis associated with a local face.
  static const int axis[6]  = { 0, 1, 2, 0, 1, 2 };

  // Gives the location of sending face on the receiver.
  static const float dir[6] = { 1, 1, 1, -1, -1, -1 };

  // Temporary store for local particle injectors.
  // FIXME: Ugly static usage.
  static particle_injector_t * RESTRICT ALIGNED(16) ci = NULL;

  static int max_ci = 0;

  #ifdef VPIC_PARTICLE_ANNOTATION
  // Same static buffer for the transport of annotations in the third run.
  // We need 6 of them for the 6 faces of the local domain through which
  // communication can happen. The sizing is for the maximum of all six.
  // cab: Communcation of Annotations Buffer
  static char * RESTRICT ALIGNED(16) cab[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
  static int max_cab = 0; // Size of the buffer (in bytes)
  int max_cas = 0; //Maximum size that might be needed for a particle (in byte)

  int sp_oldnp[MAX_SP]; // Number of particles we held per species before new particles came in.
  #endif


  int n_send[6], n_recv[6], n_ci;

  species_t* sp;

  int face;

  // Check input args.

  if ( ! sp_list )
  {
    return; // Nothing to do if no species.
  }

  if ( ! fa                ||
       ! aa                ||
       sp_list->g != aa->g ||
       fa->g      != aa->g )
  {
    ERROR( ( "Bad args." ) );
  }

  // Unpack the particle boundary conditions.

  particle_bc_func_t pbc_interact[MAX_PBC];

  void * pbc_params[MAX_PBC];

  const int nb = num_particle_bc( pbc_list );

  if ( nb > MAX_PBC )
  {
    ERROR( ( "Update this to support more particle boundary conditions." ) );
  }

  for( particle_bc_t * pbc = pbc_list; pbc; pbc = pbc->next )
  {
    pbc_interact[ -pbc->id - 3 ] = pbc->interact;
    pbc_params  [ -pbc->id - 3 ] = pbc->params;
  }

  // Unpack fields.

  field_t * RESTRICT ALIGNED(128) f = fa->f;
  grid_t  * RESTRICT              g = fa->g;

  // Unpack accumulator.

  accumulator_t * RESTRICT ALIGNED(128) a0 = aa->a;

  // Unpack the grid.

  const int64_t * RESTRICT ALIGNED(128) neighbor = g->neighbor;
  /**/  mp_t    * RESTRICT              mp       = g->mp;

  const int64_t rangel = g->rangel;
  const int64_t rangeh = g->rangeh;
  const int64_t rangem = g->range[world_size];

  /*const*/ int bc[6], shared[6];
  /*const*/ int64_t range[6];

  // debug print
  /*
  for(int r = 0; r < nproc(); r++) {
    if(r == rank()) {
      LIST_FOR_EACH( sp, sp_list ) {
        #ifdef VPIC_PARTICLE_ANNOTATION
          int sp_annotation = sp->has_annotation;
        #else
          int sp_annotation = 0;
        #endif
        printf("<%d> %d particles, %d movers, %d annotation in species %s\n", rank(), sp->np, sp->nm, sp_annotation, sp->name);
      }
    }
    mp_barrier();
  }
  */

  #ifdef VPIC_PARTICLE_ANNOTATION
  // How large might a particle be (with species ID, all annotation and possibily an ID)
  max_cas = 0;
  LIST_FOR_EACH( sp, sp_list ) {
    int cas = 0;
    // Space for the species ID
    cas += sizeof(sp->id);
    // Space for the annotations
    if(sp->has_annotation) {
       cas += sp->has_annotation*sizeof(float);
    }
    // Space for the particle ID
    #ifdef VPIC_GLOBAL_PARTICLE_ID
      cas += sizeof(sp->p_id[0]);
    #endif
    if(cas > max_cas) {
      max_cas = cas;
      /*
      if(rank() == 0) {
        printf("sp->id of size %d", sizeof(sp->id));
        if(sp->has_annotation) {
          printf(" + %d annotations of size %d", sp->has_annotation, sizeof(float));
        }
        #ifdef VPIC_GLOBAL_PARTICLE_ID
          printf(" and particle ID of size %d", sizeof(sp->p_id[0]));
        #endif
        printf(" need %d bytes\n", cas);
      }
      */
    }
  }
  // printf("<%d> particle annotations (and optional ID) need up to %d bytes per particle\n", rank(), max_cas);

  // How many max_cas sized particles might we have to communicate
  int n_cab = 0;
  LIST_FOR_EACH( sp, sp_list ) {
    n_cab += sp->nm;
  }

  // make sure the buffer we have is large enough
  // FIXME: should we shrink that ever?
  if ( max_cab < max_cas * n_cab ) {
    // printf("<%d> Reallocating to get space for %d bytes in annotations\n", rank(), max_cas*n_cab);
    for( face = 0; face < 6; face++ ) {
      char * new_cab = cab[face];
      FREE_ALIGNED( new_cab );

      MALLOC_ALIGNED( new_cab, max_cas * n_cab, 16 );

      cab[face] = new_cab;
    }
    max_cab = max_cas * n_cab;
  } else {
    for( face = 0; face < 6; face++ ) {
      memset(cab[face], 0, max_cas*n_cab);
    }
  }
  #endif

  for( face = 0; face < 6; face++ )
  {
    bc    [ face ] = g->bc[ f2b[ face ] ];

    shared[ face ] = ( bc[ face ] >= 0          ) &&
                     ( bc[ face ] <  world_size ) &&
                     ( bc[ face ] != world_rank );

    if ( shared[ face ] )
    {
      range[ face ] = g->range[ bc[ face ] ];
    }
  }

  // Begin receiving the particle counts.

  for( face = 0; face < 6; face++ )
  {
    if ( shared[ face ] )
    {
      mp_size_recv_buffer( mp,
                           f2b[ face ],
                           sizeof( int ) );

      mp_begin_recv( mp,
                     f2b[ face ],
                     sizeof( int ),
                     bc[ face ],
                     f2rb[ face ] );
    }
  }

  // Load the particle send and local injection buffers.

  do
  {
    particle_injector_t * RESTRICT ALIGNED(16) pi_send[6];

    // Presize the send and injection buffers.
    //
    // Each buffer is large enough to hold one injector corresponding
    // to every mover in use (worst case, but plausible scenario in
    // beam simulations, is one buffer gets all the movers).
    //
    // FIXME: We could be several times more efficient in our particle
    // injector buffer sizing here.  Namely, we could create one local
    // injector buffer of nm in size.  All injection for all boundaries
    // would be done here.  The local buffer would then be counted to
    // determine the size of each send buffer.  The local buffer would
    // then move all injectors into the appropriate send buffers, leaving
    // only the local injectors.  This would require some extra data
    // motion though, but would give a more robust implementation against
    // variations in MP implementation.
    //
    // FIXME: This presizing assumes that custom boundary conditions
    // inject at most one particle per incident particle.  Currently,
    // the invocation of pbc_interact[*] insures that assumption will
    // be satisfied, if the handlers conform that it.  We should be
    // more flexible though in the future, especially given the above
    // overalloc.

    int nm = 0;

    LIST_FOR_EACH( sp, sp_list ) nm += sp->nm;

    for( face = 0; face < 6; face++ )
    {
      if ( shared[ face ] )
      {
        mp_size_send_buffer( mp,
                             f2b[ face ],
                             16 + nm * sizeof( particle_injector_t ) );

        pi_send[ face ] = (particle_injector_t *) ( ( (char *) mp_send_buffer( mp,
                                                                               f2b[ face ] )
                                                    ) + 16 );

        n_send[ face ] = 0;
      }
    }

    if ( max_ci < nm )
    {
      particle_injector_t * new_ci = ci;

      FREE_ALIGNED( new_ci );

      MALLOC_ALIGNED( new_ci, nm, 16 );

      ci     = new_ci;
      max_ci = nm;
    }

    n_ci = 0;

    // For each species, load the movers.

    LIST_FOR_EACH( sp, sp_list )
    {
      const float   sp_q  = sp->q;
      const int32_t sp_id = sp->id;

      particle_t * RESTRICT ALIGNED(128) p0 = sp->p;
      #ifdef VPIC_GLOBAL_PARTICLE_ID
      int sp_has_ids = sp->has_ids;
      size_t* RESTRICT ALIGNED(128) p_id = sp->p_id; // May be NULL if sp_has_ids if false
      #endif
      #ifdef VPIC_PARTICLE_ANNOTATION
      int sp_has_annotation = sp->has_annotation;
      float* RESTRICT ALIGNED(128) sp_p_annotation = sp->p_annotation;
      #endif

      int np = sp->np;

      particle_mover_t * RESTRICT ALIGNED(16)  pm = sp->pm + sp->nm - 1;
      nm = sp->nm;

      particle_injector_t * RESTRICT ALIGNED(16) pi;
      int i, voxel;
      int64_t nn;

      // Note that particle movers for each species are processed in
      // reverse order.  This allows us to backfill holes in the
      // particle list created by boundary conditions and/or
      // communication.  This assumes particles on the mover list are
      // monotonically increasing.  That is: pm[n].i > pm[n-1].i for
      // n=1...nm-1.  advance_p and inject_particle create movers with
      // property if all aged particle injection occurs after
      // advance_p and before this.

      for( ; nm; pm--, nm-- )
      {
        i       = pm->i;
        voxel   = p0[i].i;
        face    = voxel & 7;
        voxel >>= 3;
        p0[i].i = voxel;
        nn      = neighbor[ 6 * voxel + face ];

        // Absorb.

        if ( nn == absorb_particles )
        {
          // Ideally, we would batch all rhob accumulations together
          // for efficiency.
          accumulate_rhob( f, p0 + i, g, sp_q );

          goto backfill;
        }

        // Send to a neighboring node.

        if ( ( ( nn >= 0      ) & ( nn <  rangel ) ) |
             ( ( nn >  rangeh ) & ( nn <= rangem ) ) )
        {
          pi = &pi_send[ face ] [ n_send[ face ] ];

          #ifdef V4_ACCELERATION

          copy_4x1( &pi->dx,    &p0[i].dx  );
          copy_4x1( &pi->ux,    &p0[i].ux  );
          copy_4x1( &pi->dispx, &pm->dispx );

          #else

          pi->dx    = p0[i].dx;
          pi->dy    = p0[i].dy;
          pi->dz    = p0[i].dz;

          pi->ux    = p0[i].ux;
          pi->uy    = p0[i].uy;
          pi->uz    = p0[i].uz;
          pi->w     = p0[i].w;

          pi->dispx = pm->dispx;
          pi->dispy = pm->dispy;
          pi->dispz = pm->dispz;

          #endif

          #ifdef VPIC_GLOBAL_PARTICLE_ID
          // Send global id too
          if(sp_has_ids) {
            pi->global_particle_id = p_id[i];
          }
          #endif

          #ifdef VPIC_PARTICLE_ANNOTATION
          float* RESTRICT ALIGNED(128) p_annotation = (float*) &(cab[face][n_send[face] * max_cas]);

          //printf("<%d> %dth particle of species %s has %d annotation that need to go into the buffer starting at %p\n", rank(), i, sp->name, sp_has_annotation, p_annotation);

          // Store species ID
          {
            // Oh god scary...
            int* p_annotation_int = (int*) p_annotation;
            *p_annotation_int = sp_id;
            // printf("<%d> %p set to %d\n", rank(), p_annotation_int, *p_annotation_int);
            p_annotation_int++;
            p_annotation = (float*) p_annotation_int;
          }

          // Store particle ID if needed
          #ifdef VPIC_GLOBAL_PARTICLE_ID
          // Also somewhat scary
          size_t* p_annotation_sizet = (size_t*) p_annotation;
          if(sp_has_ids) {
            *p_annotation_sizet = p_id[i];
          } else {
            *p_annotation_sizet = 0;
          }
          // printf("<%d> %p set to %dl\n", rank(), p_annotation_sizet, *p_annotation_sizet);
          p_annotation_sizet++;
          p_annotation = (float*) p_annotation_sizet;
          #endif

          if(sp_has_annotation) {
            // Store annotations
            for(int a=0; a<sp_has_annotation; a++) {
              *p_annotation = sp_p_annotation[i*sp_has_annotation+a];
              // printf("<%d> %p set to %f\n", rank(), p_annotation, *p_annotation);
              p_annotation++;
            }
          }
          #endif

          ( &pi->dx )[ axis[ face ] ] = dir[ face ];
          pi->i                       = nn - range[ face ];
          pi->sp_id                   = sp_id;

          n_send[ face ]++;

          goto backfill;
        }

        // User-defined handling.

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

        if ( ( nn >= 0  ) &
             ( nn <  nb ) )
        {
          n_ci += pbc_interact[ nn ]( pbc_params[ nn ],
                                      sp,
                                      p0 + i,
                                      pm,
                                      ci + n_ci,
                                      1,
                                      face );

          goto backfill;
        }

        // Uh-oh: We fell through.

        WARNING( ( "Unknown boundary interaction ... dropping particle "
                   "(species=%s)",
                   sp->name ) );

      backfill:

        np--;

        #if defined(V8_ACCELERATION)

        copy_8x1( &p0[i].dx, &p0[np].dx );

        #elif defined(V4_ACCELERATION)

        copy_4x1( &p0[i].dx, &p0[np].dx );
        copy_4x1( &p0[i].ux, &p0[np].ux );

        #else

        p0[i] = p0[np];

        #endif

        #ifdef VPIC_GLOBAL_PARTICLE_ID
        if(sp_has_ids) {
          p_id[i] = p_id[np];    /* keep p_id[] in sync with p[] */
        }
        #endif

        #ifdef VPIC_PARTICLE_ANNOTATION
        if(sp_has_annotation) {
          for(int a=0; a<sp_has_annotation; a++) {
            sp_p_annotation[i*sp_has_annotation+a] = sp_p_annotation[np*sp_has_annotation+a];
          }
        }
        #endif
      }

      sp->np = np;
      sp->nm = 0;
    }

  } while(0);

  #ifdef VPIC_PARTICLE_ANNOTATION
  // Figure out how many particles we currently have (per species) now that we
  // have moved all leaving particles to the buffer and before new particles
  // come in
  LIST_FOR_EACH( sp, sp_list ) {
    const int sp_id = sp->id;
    if((0 <= sp_id) && (sp_id < MAX_SP)) {
      sp_oldnp[sp_id] = sp->np;
    } else {
      ERROR( ("Invalid sp->id") );
    }
  }

  #endif

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

  for( face = 0; face < 6; face++ )
  {
    if ( shared[ face ] )
    {
      *( (int *) mp_send_buffer( mp,
                                 f2b[ face ] ) ) = n_send[ face ];

      mp_begin_send( mp,
                     f2b[ face ],
                     sizeof( int ),
                     bc[ face ],
                     f2b[ face ] );
    }
  }

  for( face = 0; face < 6; face++ )
  {
    if ( shared[ face ] )
    {
      mp_end_recv( mp,
                   f2b[ face ] );

      n_recv[ face ] = *( (int *) mp_recv_buffer( mp,
                                                  f2b[ face ] ) );

      mp_size_recv_buffer( mp,
                           f2b[ face ],
                           16 + n_recv[ face ] * sizeof( particle_injector_t ) );

      mp_begin_recv( mp,
                     f2b[ face ],
                     16 + n_recv[ face ] * sizeof( particle_injector_t ),
                     bc[ face ],
                     f2rb[ face ] );
    }
  }

  for( face = 0; face < 6; face++ )
  {
    if ( shared[ face ] )
    {
      mp_end_send( mp,
                   f2b[ face ] );

      // FIXME: ASSUMES MP WON'T MUCK WITH REST OF SEND BUFFER. IF WE
      // DID MORE EFFICIENT MOVER ALLOCATION ABOVE, THIS WOULD BE
      // ROBUSTED AGAINST MP IMPLEMENTATION VAGARIES.

      mp_begin_send( mp,
                     f2b[ face ],
                     16 + n_send[ face ] * sizeof( particle_injector_t ),
                     bc[ face ],
                     f2b[ face ] );
    }
  }

  #ifndef DISABLE_DYNAMIC_RESIZING
  // Resize particle storage to accomodate worst case inject.

  do
  {
    int n, nm;

    // Resize each species's particle and mover storage to be large
    // enough to guarantee successful injection.  If we broke down
    // the n_recv[face] by species before sending it, we could be
    // tighter on memory footprint here.

    int max_inj = n_ci;

    for( face = 0; face < 6; face++ )
    {
      if ( shared[ face ] )
      {
        max_inj += n_recv[ face ];
      }
    }

    LIST_FOR_EACH( sp, sp_list )
    {
      particle_mover_t * new_pm;
      particle_t       * new_p;
      #ifdef VPIC_GLOBAL_PARTICLE_ID
      int sp_has_ids   = sp->has_ids;
      size_t           * new_p_id;
      #endif
      #ifdef VPIC_PARTICLE_ANNOTATION
      int sp_has_annotation = sp->has_annotation;
      float                 * new_p_annotation;
      #endif

      n = sp->np + max_inj;

      if ( n > sp->max_np )
      {
        n += 0.3125 * n; // Increase by 31.25% (~<"silver
        /**/             // ratio") to minimize resizes (max
        /**/             // rate that avoids excessive heap
        /**/             // fragmentation)

        // WARNING( ( "Resizing local %s particle storage from %i to %i",
        //            sp->name,
        //            sp->max_np,
        //            n ) );

        MALLOC_ALIGNED( new_p, n, 128 );

        COPY( new_p, sp->p, sp->np );

        FREE_ALIGNED( sp->p );

        #ifdef VPIC_GLOBAL_PARTICLE_ID
        /* changes made to p[] must be mirrored in p_id[] */
        if(sp_has_ids) {
          MALLOC_ALIGNED( new_p_id, n, 128 );

          COPY( new_p_id, sp->p_id, sp->np );

          FREE_ALIGNED( sp->p_id );

          sp->p_id   = new_p_id;
        }
        #endif

        #ifdef VPIC_PARTICLE_ANNOTATION
        /* changes made to p[] must also be mirrored in p_annotation[] */
        if(sp_has_annotation) {
          MALLOC_ALIGNED( new_p_annotation, n * sp_has_annotation, 128 );

          COPY( new_p_annotation, sp->p_annotation, sp->np * sp_has_annotation);

          FREE_ALIGNED( sp->p_annotation );

          sp->p_annotation = new_p_annotation;
        }
        #endif

        sp->p      = new_p;
        sp->max_np = n;
      }

      else if( sp->max_np > MIN_NP          &&
               n          < sp->max_np >> 1 )
      {
        n += 0.125 * n; // Overallocate by less since this rank is decreasing

        if ( n < MIN_NP )
        {
          n = MIN_NP;
        }

        // WARNING( ( "Resizing (shrinking) local %s particle storage from "
        //            "%i to %i",
        //            sp->name,
        //            sp->max_np,
        //            n ) );

        MALLOC_ALIGNED( new_p, n, 128 );

        COPY( new_p, sp->p, sp->np );

        FREE_ALIGNED( sp->p );

        #ifdef VPIC_GLOBAL_PARTICLE_ID
        /* changes made to p[] must be mirrored in p_id[] */
        if(sp_has_ids) {
          MALLOC_ALIGNED( new_p_id, n, 128 );

          COPY( new_p_id, sp->p_id, sp->np );

          FREE_ALIGNED( sp->p_id );

          sp->p_id   = new_p_id;
        }
        #endif

        #ifdef VPIC_PARTICLE_ANNOTATION
        /* changes made to p[] must also be mirrored in p_annotation[] */
        if(sp_has_annotation) {
          MALLOC_ALIGNED( new_p_annotation, n * sp_has_annotation, 128 );

          COPY( new_p_annotation, sp->p_annotation, sp->np * sp_has_annotation);

          FREE_ALIGNED( sp->p_annotation );

          sp->p_annotation = new_p_annotation;
        }
        #endif

        sp->p      = new_p;
        sp->max_np = n;
      }

      // Feasibly, a vacuum-filled rank may receive a shock and need more movers
      // than available from MIN_NP.

      nm = sp->nm + max_inj;

      if ( nm > sp->max_nm )
      {
        nm += 0.3125 * nm; // See note above

        // WARNING( ( "This happened.  Resizing local %s mover storage from "
        //            "%i to %i based on not enough movers",
        //            sp->name,
        //            sp->max_nm,
        //            nm ) );

        MALLOC_ALIGNED( new_pm, nm, 128 );

        COPY( new_pm, sp->pm, sp->nm );

        FREE_ALIGNED( sp->pm );

        sp->pm     = new_pm;
        sp->max_nm = nm;
      }
    }
  } while(0);
  #endif

  do
  {
    // Unpack the species list for random acesss.

    particle_t       * RESTRICT ALIGNED(32) sp_p [ MAX_SP ];
    particle_mover_t * RESTRICT ALIGNED(32) sp_pm[ MAX_SP ];
    #ifdef VPIC_GLOBAL_PARTICLE_ID
    int sp_has_ids[ MAX_SP ];
    size_t           * RESTRICT ALIGNED(32) sp_p_id[ MAX_SP ];
    #endif
    #ifdef VPIC_PARTICLE_ANNOTATION
    int sp_has_annotation[ MAX_SP ];
    size_t           * RESTRICT ALIGNED(32) sp_p_annotation[ MAX_SP ];
    #endif

    float sp_q [ MAX_SP ];
    int   sp_np[ MAX_SP ];
    int   sp_nm[ MAX_SP ];

    #ifdef DISABLE_DYNAMIC_RESIZING
    int sp_max_np[64], n_dropped_particles[64];
    int sp_max_nm[64], n_dropped_movers   [64];
    #endif

    if ( num_species( sp_list ) > MAX_SP )
    {
      ERROR( ( "Update this to support more species." ) );
    }

    LIST_FOR_EACH( sp, sp_list )
    {
      sp_p[ sp->id ] = sp->p;
      sp_pm[ sp->id ] = sp->pm;
      sp_q [ sp->id ] = sp->q;
      sp_np[ sp->id ] = sp->np;
      sp_nm[ sp->id ] = sp->nm;

      #ifdef VPIC_GLOBAL_PARTICLE_ID
      sp_has_ids[ sp->id ] = sp->has_ids;
      if(sp_has_ids[sp->id]) {
        sp_p_id[ sp->id ] = sp->p_id;
      }
      #endif

      // The values for VPIC_PARTICLE_ANNOTATION are communicated later

      #ifdef DISABLE_DYNAMIC_RESIZING
      sp_max_np[ sp->id ] = sp->max_np;
      sp_max_nm[ sp->id ] = sp->max_nm;

      n_dropped_particles[ sp->id ] = 0;
      n_dropped_movers   [ sp->id ] = 0;
      #endif
    }

    // Inject particles.  We do custom local injection first to
    // increase message overlap opportunities.

    face = 5;

    do
    {
      particle_t                * RESTRICT ALIGNED(32) p;
      particle_mover_t          * RESTRICT ALIGNED(16) pm;
      const particle_injector_t * RESTRICT ALIGNED(16) pi;
      #ifdef VPIC_GLOBAL_PARTICLE_ID
      size_t                    * RESTRICT ALIGNED(32) p_id;
      #endif

      int np, nm, n, id;

      face++;

      if ( face == 7 )
      {
        face = 0;
      }

      if ( face == 6 )
      {
        pi = ci;
        n  = n_ci;
      }

      else if ( shared[ face ] )
      {
        mp_end_recv( mp,
                     f2b[ face ] );

        pi = (const particle_injector_t *)
             ( ( (char *) mp_recv_buffer( mp,
                                          f2b[ face ] ) ) + 16 );

        n  = n_recv[ face ];
      }

      else
      {
        continue;
      }

      // Reverse order injection is done to reduce thrashing of the
      // particle list. Particles are removed in reverse order so the
      // overall impact of removal + injection is to keep injected
      // particles in order.
      //
      // WARNING: THIS TRUSTS THAT THE INJECTORS, INCLUDING THOSE
      // RECEIVED FROM OTHER NODES, HAVE VALID PARTICLE IDS.

      pi += n - 1;

      for( ; n; pi--, n-- )
      {
        id = pi->sp_id;

        p = sp_p[id];
        np = sp_np[id];
        #ifdef VPIC_GLOBAL_PARTICLE_ID
        if(sp_has_ids[id]) {
          p_id = sp_p_id[id];
        }
        #endif

        pm = sp_pm[id];
        nm = sp_nm[id];

        #ifdef DISABLE_DYNAMIC_RESIZING
        if ( np >= sp_max_np[ id ] )
        {
          n_dropped_particles[ id ]++;

          continue;
        }
        #endif

        #ifdef V4_ACCELERATION

        copy_4x1( &p[np].dx, &pi->dx );
        copy_4x1( &p[np].ux, &pi->ux );

        #else

        p[np].dx = pi->dx;
        p[np].dy = pi->dy;
        p[np].dz = pi->dz;
        p[np].i  = pi->i;

        p[np].ux = pi->ux;
        p[np].uy = pi->uy;
        p[np].uz = pi->uz;
        p[np].w  = pi->w;

        #endif

        #ifdef VPIC_GLOBAL_PARTICLE_ID
        if(sp_has_ids[id]) {
          p_id[np] = pi->global_particle_id;

          // std::cout << "Recving particle with global_id " << pi->global_particle_id << " on rank " << _world_rank << std::endl;
        }
        #endif

        sp_np[id] = np + 1;

        #ifdef DISABLE_DYNAMIC_RESIZING
        if ( nm >= sp_max_nm[ id ] )
        {
          n_dropped_movers[ id ]++;

          continue;
        }
        #endif

        #ifdef V4_ACCELERATION

        copy_4x1( &pm[nm].dispx, &pi->dispx );

        pm[nm].i = np;

        #else

        pm[nm].dispx = pi->dispx;
        pm[nm].dispy = pi->dispy;
        pm[nm].dispz = pi->dispz;
        pm[nm].i     = np;

        #endif

        sp_nm[id] = nm + move_p( p, pm + nm, a0, g, sp_q[id] );
      }
    } while( face != 5 );

    LIST_FOR_EACH( sp, sp_list )
    {
      #ifdef DISABLE_DYNAMIC_RESIZING
      if ( n_dropped_particles[ sp->id ] )
      {
        WARNING( ( "Dropped %i particles from species \"%s\".  Use a larger "
                   "local particle allocation in your simulation setup for "
                   "this species on this node.",
                   n_dropped_particles[ sp->id ],
                   sp->name ) );
      }

      if ( n_dropped_movers[ sp->id ] )
      {
        WARNING( ( "%i particles were not completed moved to their final "
                   "location this timestep for species \"%s\".  Use a larger "
                   "local particle mover buffer in your simulation setup "
                   "for this species on this node.",
                   n_dropped_movers[ sp->id ],
                   sp->name ) );
      }
      #endif

      sp->np = sp_np[ sp->id ];
      sp->nm = sp_nm[ sp->id ];
    }

  } while(0);

  for( face = 0; face < 6; face++ )
  {
    if ( shared[ face ] )
    {
      mp_end_send( mp,
                   f2b[ face ] );
    }
  }

  #ifdef VPIC_PARTICLE_ANNOTATION
  /* if(rank() == 0) {
    // Dump buffers for debugging
    for(face = 0; face < 6; face++) {
      if ( shared[ face ] ) {
       // printf("<%d> n_send[%d] = %d\n", rank(), face, n_send[face]);
       // printf("<%d> n_recv[%d] = %d\n", rank(), face, n_recv[face]);
       printf("<%d> cab[%d] = [", rank(), face);

       for(int i = 0; i<n_send[face]; i++) {
         // Begininng of particle
         char* print_cab = &(cab[face][i*max_cas]);
         int j = 0;

         // Print species ID
         int* print_cab_int = (int*) print_cab;
         printf("%d ", *print_cab_int);
         print_cab_int++;
         j += sizeof(int);

         // Print ID if applicable
         #ifdef VPIC_GLOBAL_PARTICLE_ID
         size_t* print_cab_sizet = (size_t*) (print_cab +j);
         printf("%ld", *print_cab_sizet);
         print_cab_sizet++;
         j += sizeof(size_t);
         #endif

         //Print annotation
         for(; j<max_cas; j+=sizeof(float)) {
           float* print_cab_float = (float*) (print_cab + j);
           printf(" %f", *print_cab_float);
         }

         printf(", ");
       }
       printf("]\n");
      }
    }
  } */

  // Ensure sufficent buffer size and prepost recv
  for( face = 0; face < 6; face++ ) {
    if ( shared[ face ] ) {
      if(n_recv[face] > 0) {
        mp_size_recv_buffer( mp, f2b[ face ], max_cas * n_recv[face] );

        mp_begin_recv( mp, f2b[ face ], max_cas * n_recv[face], bc[ face ], f2rb[ face ] );
      }
      if(n_send[face] > 0) {
        mp_size_send_buffer( mp, f2b[ face ], max_cas * n_send[face] );
      }
    }
  }

  // Fill and send buffers
  for( face = 0; face < 6; face++ ) {
    if ( shared[ face ] ) {
      if(n_send[face] > 0) {
        float* send_buf = (float *) mp_send_buffer( mp, f2b[ face ] );

        memcpy(send_buf, cab[face], max_cas*n_send[face]);

        mp_begin_send( mp, f2b[ face ], max_cas*n_send[face], bc[ face ], f2b[ face ] );
      }
    }
  }

  // Wait for recv and handle buffer
  for( face = 0; face < 6; face++ ) {
    if ( shared[ face ] ) {
      if(n_recv[face] > 0) {
        mp_end_recv( mp, f2b[face] );

        /* if(rank() == 0) {
          // Print for debugging purposes
          printf("<%d> recv_buffer[%d] = [", rank(), face);
          for(int i = 0; i<n_recv[face]; i++) {
            char* recv_buffer = (char*) mp_recv_buffer(mp,f2b[face]);
            // Begininng of particle
            char* print_cab = &(recv_buffer[i*max_cas]);
            int j = 0;

            //Print species ID
            int* print_cab_int = (int*) print_cab;
            printf("%d ", *print_cab_int);
            print_cab_int++;
            j += sizeof(int);

            // Print ID if applicable
            #ifdef VPIC_GLOBAL_PARTICLE_ID
            size_t* print_cab_sizet = (size_t*) (print_cab +j);
            printf("%ld", *print_cab_sizet);
            print_cab_sizet++;
            j += sizeof(size_t);
            #endif

            //Print annotation
            for(; j<max_cas; j+=sizeof(float)) {
              float* print_cab_float = (float*) (print_cab + j);
              printf(" %f", *print_cab_float);
            }

            printf(", ");
          }
          printf("]\n");
        } */

        // Unpack buffer and store attributes
        for(int i = n_recv[face]-1; i>=0; i--) {
          char* recv_buffer = ((char*) mp_recv_buffer(mp,f2b[face])) + i*max_cas ;
          // Get species ID
          const int sp_id = *( (int*) recv_buffer );
          recv_buffer += sizeof(int);
          // Get particle ID if applicable
          #ifdef VPIC_GLOBAL_PARTICLE_ID
          const size_t p_id = *( (size_t*) recv_buffer );
          recv_buffer += sizeof(size_t);
          #endif
          // The rest of the buffer must be annotations
          const float* p_annotation = (float*) recv_buffer;

          species_t * sp = find_species_id(sp_id,sp_list); // We should probably do that once for all possible MAX_SP ids
          if(!sp) {
             printf("Can not find species for sp_id = %d", sp_id);
             printf("Particle %d of %d on face %d\n", i, n_recv[face], face);
             ERROR();
          }
          #ifdef VPIC_GLOBAL_PARTICLE_ID
          if(sp->has_ids) {
            if(sp->p_id[sp_oldnp[sp_id]] != p_id) {
              printf("We are at particle %d of species %d which has ID %ld, not %ld as received\n", sp_oldnp[sp_id], sp_id, sp->p_id[sp_oldnp[sp_id]], p_id);
              ERROR( ("ID missmatch") );
            } else {
              // printf("We are at particle %d of species %d which has ID %ld, matching the received %ld\n", sp_oldnp[sp_id], sp_id, sp->p_id[sp_oldnp[sp_id]], p_id);
            }
          }
          #endif

          if(sp->has_annotation) {
            for(int a=0; a<sp->has_annotation; a++) {
              sp->p_annotation[sp_oldnp[sp_id]*sp->has_annotation + a] = p_annotation[a];
            }
          }

          // The next particle that comes in will be at the next index
          sp_oldnp[sp_id]++;
        }
      }
    }
  }

  // Wait for send to end
  for( face = 0; face < 6; face++ ) {
    if ( shared[ face ] ) {
      if(n_send[face] > 0) {
        mp_end_send( mp, f2b[face] );
      }
    }
  }
  #endif
}
