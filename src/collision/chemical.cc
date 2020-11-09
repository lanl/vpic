#define IN_collision

#include "chemical.h"

/* Private interface *********************************************************/

void
apply_chemical_collision_model( chemical_collision_model_t * cm )
{

  species_t * sp;
  species_t ** reactants = cm->reactants;
  grid_t    *  g = cm->reactants[0]->g;

  const int N = cm->n_reactants;
  int n;

  if ( cm->interval < 1 || ( g->step % cm->interval ) )
  {
    return;
  }

  FOR_REACTANTS if( sp->last_sorted != g->step ) sort_p( sp );

  // Conditionally execute this when more abstractions are available.
  apply_chemical_collision_model_pipeline( cm );

}

void
checkpt_chemical_collision_model( const collision_op_t * cop ) {
  const chemical_collision_model_t * cm =
    (const chemical_collision_model_t *)cop->params;
  int i;
  CHECKPT( cm, 1 );
  CHECKPT_STR( cm->name );
  CHECKPT_SYM( cm->rate_constant );
  CHECKPT_SYM( cm->collision );
  CHECKPT_PTR( cm->params );
  CHECKPT_PTR( cm->rp );
  CHECKPT_PTR( cm->fa );
  CHECKPT( cm->consumable, cm->n_reactants );
  for(i=0 ; i < cm->n_products ; ++i ) CHECKPT_PTR( cm->products[i] );
  for(i=0 ; i < cm->n_reactants ; ++i ) CHECKPT_PTR( cm->reactants[i] );

  checkpt_collision_op_internal( cop );
}

collision_op_t *
restore_chemical_collision_model( void ) {
  chemical_collision_model_t * cm;
  int i;
  RESTORE( cm );
  MALLOC( cm->reactants,         cm->n_reactants );
  MALLOC( cm->products,          cm->n_products  );
  MALLOC( cm->n_produced,                          MAX_PIPELINE );
  MALLOC( cm->n_modified,        cm->n_reactants * MAX_PIPELINE );
  MALLOC( cm->reactant_movers,   cm->n_reactants * MAX_PIPELINE );
  MALLOC( cm->product_particles, cm->n_products  * MAX_PIPELINE );
  RESTORE_STR( cm->name );
  RESTORE_SYM( cm->rate_constant );
  RESTORE_SYM( cm->collision );
  RESTORE_PTR( cm->params );
  RESTORE_PTR( cm->rp );
  RESTORE_PTR( cm->fa );
  RESTORE( cm->consumable );
  for(i=0 ; i < cm->n_products ; ++i ) RESTORE_PTR( cm->products[i] );
  for(i=0 ; i < cm->n_reactants ; ++i ) RESTORE_PTR( cm->reactants[i] );
  return restore_collision_op_internal( cm );
}

void
delete_chemical_collision_model( collision_op_t * cop ) {
  chemical_collision_model_t * cm = (chemical_collision_model_t *)cop->params;
  FREE( cm->name );
  FREE( cm->reactants ); FREE( cm->reactant_movers );   FREE( cm->n_modified );
  FREE( cm->products  ); FREE( cm->product_particles ); FREE( cm->n_produced );
  FREE( cm->consumable );
  FREE( cm );
  delete_collision_op_internal( cop );
}

/* Public interface **********************************************************/

collision_op_t *
chemical_collision_model( const char               * RESTRICT name,
                          /**/  void               * RESTRICT params,
                          /**/  species_t         ** RESTRICT reactants,
                          const int                * RESTRICT consumable,
                          const int                           n_reactants,
                          /**/  species_t         ** RESTRICT products,
                          const int                           n_products,
                          /**/  chemical_rate_constant_func_t rate_constant,
                          /**/  chemical_collision_func_t     collision,
                          /**/  rng_pool_t         * RESTRICT rp,
                          /**/  field_array_t      * RESTRICT fa,
                          const double                        sample,
                          const int                           dynamic_sampling,
                          const int                           interval )
{
  chemical_collision_model_t * cm;
  size_t len = name ? strlen(name) : 0;

  // if( !rate_constant || !collision || !spi || !spj || spi->g!=spj->g ||
  //     !rp || rp->n_rng<N_PIPELINE ) ERROR(( "Bad args" ));

  if( len==0 )
    ERROR(( "Cannot specify a nameless collision model" ));

  if( params && !object_id( params ) )
    ERROR(( "collision model parameters must be checkpoint registered" ));

  if( n_reactants <= 0 || n_products <= 0 )
    ERROR(("Reactants and products must be >= 1 for chemical collisions."));

  const int N = n_reactants;
  const int M = n_products;

  MALLOC( cm, 1 );
  MALLOC( cm->name, len+1 );  strcpy( cm->name, name );
  MALLOC( cm->reactants,  N );  COPY( cm->reactants,  reactants,  N);
  MALLOC( cm->consumable, N );  COPY( cm->consumable, consumable, N);
  MALLOC( cm->products,   M  ); COPY( cm->products,   products,   M);
  MALLOC( cm->n_produced,          MAX_PIPELINE );
  MALLOC( cm->n_modified,        N*MAX_PIPELINE );
  MALLOC( cm->reactant_movers,   N*MAX_PIPELINE );
  MALLOC( cm->product_particles, M*MAX_PIPELINE );

  cm->n_reactants      = N;
  cm->n_products       = M;
  cm->dynamic_sampling = dynamic_sampling;
  cm->rate_constant    = rate_constant;
  cm->collision        = collision;
  cm->params           = params;
  cm->rp               = rp;
  cm->fa               = fa;
  cm->sample           = sample;
  cm->interval         = interval;

  return new_collision_op_internal( cm,
                                    (collision_op_func_t)apply_chemical_collision_model,
                                    delete_chemical_collision_model,
                                    (checkpt_func_t)checkpt_chemical_collision_model,
                                    (restore_func_t)restore_chemical_collision_model,
                                    NULL );
}
