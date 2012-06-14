/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include <material.h>

/* Private interface *********************************************************/

void
checkpt_material( const material_t * m ) {
  CHECKPT( m, 1 );
  CHECKPT_STR( m->name );
  CHECKPT_PTR( m->next );
}

material_t *
restore_material( void ) {
  material_t * m;
  RESTORE( m );
  RESTORE_STR( m->name );
  RESTORE_PTR( m->next );
  return m;
}

void
delete_material( material_t * m ) {
  UNREGISTER_OBJECT( m );
  FREE( m->name );
  FREE( m );
}

/* Public interface **********************************************************/

int
num_material( const material_t * m_list ) {
  return m_list ? m_list->id+1 : 0;
}

void
delete_material_list( material_t * m_list ) {
  material_t * m;
  while( m_list ) {
    m = m_list;
    m_list = m_list->next;
    delete_material( m );
  }
}

material_t *
find_material_name( const char * name, 
                    material_t * m_list ) {
  material_t * m;
  if( !name ) return NULL;
  LIST_FIND_FIRST( m, m_list, strcmp( m->name, name )==0 );
  return m;
}

material_t *
find_material_id( material_id id, 
                  material_t * m_list ) {
  material_t * m;
  LIST_FIND_FIRST( m, m_list, m->id==id );
  return m;
}

material_id
get_material_id( const material_t * m ) {
  return m ? m->id : -1;
}

material_t *
append_material( material_t * m,
                 material_t ** m_list ) {
  int id;
  if( !m || !m_list ) ERROR(( "Bad args" ));
  if( m->next ) ERROR(( "Material \"%s\" already in a list", m->name ));
  if( find_material_name( m->name, *m_list ) )
    ERROR(( "There is already a material named \"%s\" in list", m->name ));
  id = num_material( *m_list );
  if( id>=max_material )
    ERROR(( "Too many materials in list to append material \"%s\"", m->name ));
  m->id   = (material_id)id;
  m->next = *m_list;
  *m_list = m;
  return m;
}

material_t *
material( const char * name,
          float epsx,   float epsy,   float epsz,
          float mux,    float muy,    float muz,
          float sigmax, float sigmay, float sigmaz,
          float zetax,  float zetay,  float zetaz ) {
  material_t *m;
  int len = name ? strlen( name ) : 0;
  if( !len ) ERROR(( "Cannot create a nameless material" ));
  MALLOC( m, 1 );
  CLEAR( m, 1 );
  MALLOC( m->name, len+1 );
  strcpy( m->name, name );
  m->epsx   = epsx,   m->epsy   = epsy,   m->epsz   = epsz;
  m->mux    = mux,    m->muy    = muy,    m->muz    = muz;
  m->sigmax = sigmax, m->sigmay = sigmay, m->sigmaz = sigmaz;
  m->zetax  = zetax,  m->zetay  = zetay,  m->zetaz  = zetaz;
  /* id, next set by append_material */
  REGISTER_OBJECT( m, checkpt_material, restore_material, NULL );
  return m;
}

