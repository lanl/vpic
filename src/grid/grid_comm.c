#include "grid.h"

// FIXME: THIS CODE IS PRETTY SILLY AT THIS POINT

// FIXME: REMOTE PORT / LOCAL PORT HOOKUPS SHOULD BE EXPLICITLY
// SPECIFIED IN THE GRID.

// FIXME: THE PORT SPECIFIED IN THE API FOR BEGIN_,END_RECV IS THE
// REMOTE _SENDING_ PORT.  HOW UGLY.  USE THE RECV PORT HERE AND
// UPDATE THE CALLERS TO NOT SUCK.

// FIXME: THE grid_t SHOULD BE THE FIRST ARG PASSED (NOT THE LAST)

// Note: Messages are tagged by the port on the sender

void
begin_recv_port( int i, int j, int k,
                 int size,
                 const grid_t * g ) {
  int port = BOUNDARY(-i,-j,-k), src = g->bc[port];
  if( src<0 || src>=world_size ) return;
  mp_size_recv_buffer( g->mp, BOUNDARY(-i,-j,-k), size );
  mp_begin_recv( g->mp, port, size, src, BOUNDARY(i,j,k) );
}

void * ALIGNED(128)
end_recv_port( int i, int j, int k,
               const grid_t * g ) {
  int port = BOUNDARY(-i,-j,-k), src = g->bc[port];
  if( src<0 || src>=world_size ) return NULL;
  mp_end_recv( g->mp, port );
  return mp_recv_buffer( g->mp, port );
}

void * ALIGNED(128)
size_send_port( int i, int j, int k,
                int size,
                const grid_t * g ) {
  int port = BOUNDARY( i, j, k), dst = g->bc[port];
  if( dst<0 || dst>=world_size ) return NULL;
  mp_size_send_buffer( g->mp, port, size );
  return mp_send_buffer( g->mp, port );
}

void
begin_send_port( int i, int j, int k,
                 int size,
                 const grid_t * g ) {
  int port = BOUNDARY( i, j, k), dst = g->bc[port];
  if( dst<0 || dst>=world_size ) return;
  mp_begin_send( g->mp, port, size, dst, port );
}

void
end_send_port( int i, int j, int k,
               const grid_t * g ) {
  int port = BOUNDARY( i, j, k), dst = g->bc[port];
  if( dst<0 || dst>=world_size ) return;
  mp_end_send( g->mp, BOUNDARY(i,j,k) );
}

