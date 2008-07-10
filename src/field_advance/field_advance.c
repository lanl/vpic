#include "field_advance.h"

field_advance_t *
new_field_advance( grid_t * g,
                   material_t * m_list,
                   field_advance_methods_t * fam ) {
  field_advance_t * fa;

  // FIXME: MORE ROBUST CHECKING OF FAM
  if( g==NULL || m_list==NULL || fam==NULL )  ERROR(( "Bad args" ));

  MALLOC( fa, 1 );

  fa->method[0] = fam[0];
  fa->f = fa->method->new_field( g );
  fa->m = fa->method->new_material_coefficients( g, m_list );
  fa->g = g;

  return fa;
}

void
delete_field_advance( field_advance_t * fa ) {
  if( fa==NULL ) return; // Do nothing request
  fa->method->delete_material_coefficients( fa->m );
  fa->method->delete_field( fa->f );
  FREE( fa );
}
