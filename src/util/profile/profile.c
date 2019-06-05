#include "profile.h"
#include "sys/time.h"

profile_internal_use_only_timer_t profile_internal_use_only[] = {
# define PROFILE_TIMER_INIT( timer ) { #timer, 0., 0., 0, 0 },
  PROFILE_TIMERS( PROFILE_TIMER_INIT )
# undef PROFILE_TIMER_INIT
  { NULL, 0., 0., 0, 0 }
};

double
update_profile( int dump ) {
  profile_internal_use_only_timer_t * p;
  double sum = 0, sum_total = 0;

  for( p=profile_internal_use_only; p->name; p++ ) {
    p->t_total += p->t;
    p->n_total += p->n;
    sum        += p->t;
    sum_total  += p->t_total;
  }

  if( dump )
  {
    #if defined(VPIC_PRINT_MORE_DIGITS)
    log_printf( "\n" // 8901234567890123456 | xxx% x.xxxe+xx x.xe+xx x.xxxe+xx | xxx% x.xxxe+xx x.xe+xx x.xxxe+xx
                "                           |      Since   Last Update         |      Since   Last Restore\n"
                "    Operation              | Pct   Time      Count      Per   | Pct   Time      Count      Per\n"
                "---------------------------+----------------------------------+----------------------------------\n" );
    #else
    log_printf( "\n" // 8901234567890123456 | xxx% x.xe+xx x.xe+xx x.xe+xx | xxx% x.xe+xx x.xe+xx x.xe+xx
                "                           |      Since Last Update       |     Since Last Restore\n"
                "    Operation              | Pct   Time    Count    Per   | Pct   Time    Count    Per\n"
                "---------------------------+------------------------------+------------------------------\n" );
    #endif

    for( p=profile_internal_use_only; p->name; p++ ) {
      if( p->n==0 && p->n_total==0 ) continue;
      #if defined(VPIC_PRINT_MORE_DIGITS)
      log_printf( "%26.26s | % 3d%% %.3e %.1e %.3e | % 3d%% %.3e %.1e %.3e\n",
      #else
      log_printf( "%26.26s | % 3d%% %.1e %.1e %.1e | % 3d%% %.1e %.1e %.1e\n",
      #endif
                  p->name,
                  (int)( 100.*p->t/sum + 0.5 ), p->t,
                  (double)p->n,
                  p->t/(DBL_EPSILON+(double)p->n ),
                  (int)( 100.*p->t_total/sum_total + 0.5 ), p->t_total,
                  (double)p->n_total,
                  p->t_total/(DBL_EPSILON+(double)p->n_total) );
    }

    log_printf( "\n" );
  }

  for( p=profile_internal_use_only; p->name; p++ ) {
    p->t = 0;
    p->n = 0;
  }
  return sum;
}

double
wallclock( void ) {
  struct timeval tv[1];
  gettimeofday( tv, NULL );
  return (double)(tv->tv_sec) + 1e-6*(double)(tv->tv_usec);
}
