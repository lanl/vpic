#define IN_collision
#include <collision_private.h>

/* Private interface *********************************************************/

void
checkpt_collision_op_internal( const collision_op_t * RESTRICT cop ) {
  CHECKPT( cop, 1 );
  CHECKPT_SYM( cop->apply );
  CHECKPT_SYM( cop->delete_cop );
  CHECKPT_PTR( cop->next );
}

collision_op_t *
restore_collision_op_internal( void * params ) {
  collision_op_t * cop;
  RESTORE( cop );
  cop->params = params;
  RESTORE_SYM( cop->apply );
  RESTORE_SYM( cop->delete_cop );
  RESTORE_PTR( cop->next );
  return cop;
}

collision_op_t *
new_collision_op_internal( void * params,
                           collision_op_func_t apply,
                           delete_collision_op_func_t delete_cop,
                           checkpt_func_t checkpt,
                           restore_func_t restore,
                           reanimate_func_t reanimate ) {
  collision_op_t * cop;
  MALLOC( cop, 1 );
  cop->params     = params;
  cop->apply      = apply;
  cop->delete_cop = delete_cop;
  cop->next       = NULL; /* Set by append_collision_op */
  REGISTER_OBJECT( cop, checkpt, restore, reanimate );
  return cop;
}

void
delete_collision_op_internal( collision_op_t * cop ) {
  UNREGISTER_OBJECT( cop );
  FREE( cop );
}

/* Public interface **********************************************************/

int
num_collision_op( const collision_op_t * RESTRICT cop_list ) {
  const collision_op_t * RESTRICT cop;
  int n = 0;
  LIST_FOR_EACH( cop, cop_list ) n++;
  return n;
}

void
apply_collision_op_list( collision_op_t * cop_list ) {
  collision_op_t * cop;
  LIST_FOR_EACH( cop, cop_list ) cop->apply( cop->params );
}

void
delete_collision_op_list( collision_op_t * cop_list ) {
  collision_op_t * cop;
  while( cop_list ) {
    cop = cop_list;
    cop_list = cop_list->next;
    cop->delete_cop( cop );
  }
}

collision_op_t *
append_collision_op( collision_op_t * cop,
                     collision_op_t ** cop_list ) {
  if( !cop || !cop_list || cop->next ) ERROR(( "Bad args" ));
  cop->next = *cop_list;
  *cop_list = cop;
  return cop;
}

