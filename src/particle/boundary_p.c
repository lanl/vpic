/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include <species.h>    /* For species_t */ 
#include <mtrand.h>     /* For mt_handle */ 

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,g->nx+1,0,g->ny+1,0,g->nz+1)]
#define CUSTOM_PBC_MIN_INJECTORS 16 
#define CUSTOM_PBC_RESIZE_FACTOR 1.1

int boundary_p( particle_mover_t * ALIGNED pm,
                int nm,
		int max_nm,
		particle_t * ALIGNED p,
		int np,
		int max_np,
		field_t * ALIGNED f,
		accumulator_t * ALIGNED a,
		const grid_t * g,
		species_t * sp,
		mt_handle rng ) {
  const int sf2b[6] = { BOUNDARY(-1, 0, 0),
                        BOUNDARY( 0,-1, 0),
                        BOUNDARY( 0, 0,-1),
                        BOUNDARY( 1, 0, 0),
                        BOUNDARY( 0, 1, 0),
                        BOUNDARY( 0, 0, 1) };
  const int rf2b[6] = { BOUNDARY( 1, 0, 0),
                        BOUNDARY( 0, 1, 0),
                        BOUNDARY( 0, 0, 1),
                        BOUNDARY(-1, 0, 0),
                        BOUNDARY( 0,-1, 0),
                        BOUNDARY( 0, 0,-1) };
  int rank, nproc, pass, face, ns[6], n, *buf;
  particle_injector_t *ps[6], *pr;
  error_code err;
  static int cpb_size=CUSTOM_PBC_MIN_INJECTORS; 
  static particle_injector_t *cmlist=NULL, *cm; 

  if( pm==NULL ) ERROR(("Bad particle mover"));
  if( nm<0     ) ERROR(("Bad number of movers"));
  if( max_nm<0 ) ERROR(("Bad max number of movers"));
  if( p==NULL  ) ERROR(("Bad particle array"));
  if( np<0     ) ERROR(("Bad number of particles"));
  if( max_np<0 ) ERROR(("Bad max number of particles"));
  if( f==NULL  ) ERROR(("Bad field"));
  if( a==NULL  ) ERROR(("Bad accumulator"));
  if( g==NULL  ) ERROR(("Bad grid"));

  rank = mp_rank(g->mp);
  nproc = mp_nproc(g->mp);

# define SHARED_REMOTELY(bound) \
  (g->bc[bound]>=0 && g->bc[bound]<nproc && g->bc[bound]!=rank)

/* 3 passes may no longer suffice with custom pbcs like reflux or secondary 
   emission in multi-dimensions. */ 

# define MAX_PASSES 6   

  for( pass=0; pass<MAX_PASSES; pass++ ) {  

    /* Create particle injectors ... the sizing is probably overkill ... it
       ensures that any given send buffer is large enough to hold injectors
       for the entire mover list */
    for( face=0; face<6; face++ ) {
      ns[face] = 0;
      if( SHARED_REMOTELY(sf2b[face]) ) {
        err = mp_size_send_buffer( sf2b[face],
                                   sizeof(int) +
                                   nm*sizeof(particle_injector_t),
                                   g->mp );
        if( err!=SUCCESS ) ERROR(("Unable to size send buffer - %s",err));
        buf = (int *)mp_send_buffer(sf2b[face],g->mp);
        ps[face] = (particle_injector_t *)(&buf[1]);
      }
    }

    /* Create space for particle injectors created from custom pbc */   
    if ( !cmlist ) {
      cpb_size=( nm<CUSTOM_PBC_MIN_INJECTORS ? CUSTOM_PBC_MIN_INJECTORS : nm ); 
      cmlist=(particle_injector_t *)malloc((size_t)cpb_size*sizeof(particle_injector_t)); 
      if ( !cmlist ) ERROR(("Could not alloate space for cpb injector array.")); 
    } else if ( nm>cpb_size ) {  /* Need to resize injector array */ 
      particle_injector_t *cmlist_tmp; 
      cpb_size=nm*CUSTOM_PBC_RESIZE_FACTOR;
      cmlist_tmp=realloc( cmlist, (size_t)cpb_size*sizeof(particle_injector_t) ); 
      if ( !cmlist_tmp ) ERROR(("Could not realloate space for cpb injector array.")); 
      cmlist=cmlist_tmp; 
    }

    cm = cmlist;

    /* Load the particle send buffers
       Note: particle mover is processed in reverse order. This allows us to
       backfill holes in the particle list created by absorption and/or 
       communication. This assumes particles on the mover list are
       monotonically increasing. That is: pm[n].i > pm[n-1].i for n=1...nm-1.
       advance_p and inject_particle create movers with property if all aged
       particle injection occurs after advance_p and before this */

    for( ; nm; nm-- ) {
      particle_t *r = p + pm[nm-1].i;
#     define TEST_FACE(FACE,cond)                                       \
      if(cond) {                                                        \
        int64_t nn = g->neighbor[ 6*r->i + FACE ];                      \
        if( nn==absorb_particles ) goto done_testing;                   \
        if( (nn>=0 && nn<g->range[rank]) ||                             \
            (nn>=g->range[rank+1] && nn<g->range[nproc] ) ) {           \
          ps[FACE]->dx    = (FACE==0 || FACE==3) ? -r->dx : r->dx;      \
          ps[FACE]->dy    = (FACE==1 || FACE==4) ? -r->dy : r->dy;      \
          ps[FACE]->dz    = (FACE==2 || FACE==5) ? -r->dz : r->dz;      \
          ps[FACE]->i      = nn - g->range[g->bc[sf2b[FACE]]];          \
          /* nn - g->range[...] is less than 2^31 / 6 */                \
          ps[FACE]->ux    = r->ux;                                      \
          ps[FACE]->uy    = r->uy;                                      \
          ps[FACE]->uz    = r->uz;                                      \
          ps[FACE]->q     = r->q;                                       \
          ps[FACE]->dispx = pm[nm-1].dispx;                             \
          ps[FACE]->dispy = pm[nm-1].dispy;                             \
          ps[FACE]->dispz = pm[nm-1].dispz;                             \
          ps[FACE]++;                                                   \
          ns[FACE]++;                                                   \
          goto done_testing;                                            \
	}                                                               \
                                                                        \
        /* Test for user-defined boundary handler.          */          \
        /* Note: Particle is destroyed after it is handled. */          \
                                                                        \
        nn = -nn - 3; /* Assumes reflecting/absorbing are -1, -2 */     \
        if ( nn>=0 && nn<g->nb ) {                                      \
          g->boundary[nn].handler( g->boundary[nn].params,              \
                                   r, &pm[nm-1], f, a, g, sp, &cm, rng, FACE ); \
          goto done_testing;                                            \
        }                                                               \
      }

      TEST_FACE(0,r->dx==-1 && r->ux<0);
      TEST_FACE(1,r->dy==-1 && r->uy<0);
      TEST_FACE(2,r->dz==-1 && r->uz<0);
      TEST_FACE(3,r->dx== 1 && r->ux>0);
      TEST_FACE(4,r->dy== 1 && r->uy>0);
      TEST_FACE(5,r->dz== 1 && r->uz>0);
#     undef TEST_FACE
      WARNING(("Unknown boundary interaction ... using absorption"));
      WARNING(("pass=%i, species=%i, rank=%i", pass, sp->id, mp_rank(g->mp)));
    done_testing:
      np -= remove_p( r, p, np, f, g );
    }

    /* Exchange particle counts */

    for( face=0; face<6; face++ )
      if( SHARED_REMOTELY(rf2b[face]) ) {
        err = mp_size_recv_buffer( rf2b[face], sizeof(int), g->mp );
        if( err!=SUCCESS ) ERROR(("Unable to size recv buffer - %s",err));
        mp_begin_recv( rf2b[face], sizeof(int),
                       g->bc[rf2b[face]], sf2b[face], g->mp );
      }
    for( face=0; face<6; face++ )
      if( SHARED_REMOTELY(sf2b[face]) ) {
        buf = (int *)mp_send_buffer( sf2b[face], g->mp );
        buf[0] = ns[face];
        mp_begin_send( sf2b[face], sizeof(int),
                       g->bc[sf2b[face]], sf2b[face], g->mp );
      }
    for( face=0; face<6; face++ )
      if( SHARED_REMOTELY(rf2b[face]) ) mp_end_recv( rf2b[face], g->mp );
    for( face=0; face<6; face++ )
      if( SHARED_REMOTELY(sf2b[face]) ) mp_end_send( sf2b[face], g->mp );
    
    /* Exchange particles */

    for( face=0; face<6; face++ )
      if( SHARED_REMOTELY(rf2b[face]) ) {
        buf = (int *)mp_recv_buffer( rf2b[face], g->mp );
        n = sizeof(int) + buf[0]*sizeof(particle_injector_t);
        err = mp_size_recv_buffer( rf2b[face], n, g->mp );
        if( err!=SUCCESS ) ERROR(("Unable to size recv buffer - %s",err));
        mp_begin_recv( rf2b[face], n,
                       g->bc[rf2b[face]], sf2b[face], g->mp );
      }
    for( face=0; face<6; face++ )
      if( SHARED_REMOTELY(sf2b[face]) )
        mp_begin_send( sf2b[face],
                       sizeof(int) + ns[face]*sizeof(particle_injector_t),
                       g->bc[sf2b[face]], sf2b[face], g->mp );
    for( face=0; face<6; face++ )
      if( SHARED_REMOTELY(rf2b[face]) ) {
        mp_end_recv( rf2b[face], g->mp );
        buf = (int *)mp_recv_buffer( rf2b[face], g->mp );
        n = buf[0];
        pr = (particle_injector_t *)&buf[1];
        for(;n;n--) {
          if( np==max_np || nm==max_nm ) break;
          nm += inject_p( p, np, pm+nm, f, a, pr, g );
          pr++;
          np++;
        }
        if( n!=0 )
          WARNING(( "Ran out of room for remote injection on proc %i "
                    "(%i movers left, np = %i, max_np = %i, nm = %i, max_nm = %i)",
                    rank, n, np, max_np, nm, max_nm ));
      }
    for( face=0; face<6; face++ )
      if( SHARED_REMOTELY(sf2b[face]) ) mp_end_send( sf2b[face], g->mp );

    /* Handle particle injectors from custom pbcs */ 
    while ( cm!=cmlist ) {
      if ( np==max_np || nm==max_nm ) break;
      nm += inject_p( p, np, pm+nm, f, a, --cm, g );
      np++; 
    } 
    if ( cm!=cmlist ) 
      WARNING(( "Ran out of room for custom boundary injection on proc %i (%i movers)", 
                rank, cm-cmlist)); 
  }

  if( nm>0 ) WARNING(("Ignoring %i unprocessed movers on rank %i", nm, rank));
  return np; /* New number of particles */
}
