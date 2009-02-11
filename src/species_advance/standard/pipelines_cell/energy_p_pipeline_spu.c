#define IN_spa
#include <spa_private.h>
#include <stdio.h>

void
_SPUEAR_energy_p_pipeline_spu ( MEM_PTR( sort_p_pipeline_args_t, 128 ) argp,
                                int pipeline_rank,
								int n_pipeline ) {
  fprintf( stdout, "In energy_p_pipeline_spu\n" ); fflush( stdout );
} // _SPUEAR_energy_p_pipeline_spu
