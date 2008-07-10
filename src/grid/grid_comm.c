#include "grid.h"

// FIXME: REMOTE PORT / LOCAL PORT HOOKUPS SHOULD BE EXPLICITLY
// SPECIFIED IN THE GRID.

void
begin_recv_port( int i,
                 int j,
                 int k,
                 int size,
                 const grid_t * g ) {
  int sbound, rbound, sender;

  sbound = BOUNDARY(i,j,k);
  rbound = BOUNDARY(-i,-j,-k);
  sender = g->bc[ rbound ];
  if( sender<0 || sender>=mp_nproc(g->mp) ) return;
  mp_size_recv_buffer( sbound, size, g->mp );
  mp_begin_recv( sbound, size, sender, sbound, g->mp );
  return;
}

void * ALIGNED(16)
size_send_port( int i,
                int j,
                int k,
                int size,
                const grid_t * g ) {
  int sbound, receiver;

  sbound = BOUNDARY(i,j,k);
  receiver = g->bc[sbound];
  if( receiver<0 || receiver>=mp_nproc(g->mp) ) return NULL;
  mp_size_send_buffer( sbound, size, g->mp );
  return (void * ALIGNED(16))mp_send_buffer( sbound, g->mp );
}

void
begin_send_port( int i,
                 int j,
                 int k,
                 int size,
                 const grid_t * g ) {
  int sbound, receiver;

  sbound = BOUNDARY(i,j,k);
  receiver = g->bc[sbound];
  if( receiver<0 || receiver>=mp_nproc(g->mp) ) return;
  mp_begin_send( sbound, size, receiver, sbound, g->mp );
}

void * ALIGNED(16)
end_recv_port( int i,
               int j,
               int k,
               const grid_t * g ) {
  int sbound, rbound, sender;

  sbound = BOUNDARY(i,j,k);
  rbound = BOUNDARY(-i,-j,-k);
  sender = g->bc[rbound];
  if( sender<0 || sender>=mp_nproc(g->mp) ) return NULL;
  mp_end_recv( sbound, g->mp );
  return (void * ALIGNED(16))mp_recv_buffer(sbound,g->mp);
}

void
end_send_port( int i,
               int j,
               int k,
               const grid_t * g ) {
  int sbound, receiver;

  sbound = BOUNDARY(i,j,k);
  receiver = g->bc[sbound];
  if( receiver<0 || receiver>=mp_nproc(g->mp) ) return;
  mp_end_send( sbound, g->mp );
}
