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

/* Note: new_material added the created species to the head of the species
 * list. Further the species ids are simply incremented from the previous
 * head of the list. The first species is numbered 0. As a result, the total
 * number of species in the species list is the id of the species at the
 * head of the list plus one. */

int num_materials( const material_t *m_list ) {
  if( m_list==NULL ) return 0;
  return m_list->id+1;
}

material_id new_material( const char *name,
			  float epsx,   float epsy,   float epsz,
			  float mux,    float muy,    float muz,
			  float sigmax, float sigmay, float sigmaz,
			  material_t **m_list ) {
  material_t *m;
  material_id id;
  int len;

  if( m_list==NULL ) ERROR(("Invalid material list."));
  /* Note: strlen does not include terminating NULL */
  len = (name==NULL) ? 0 : strlen(name);
  if( len<=0 ) ERROR(("Cannot create a nameless material."));
  if( find_material_name(name,*m_list)!=NULL )
    ERROR(("There is already a material named \"%s\".", name));
  id = num_materials(*m_list);
  if( id==max_num_materials )
    ERROR(("Cannot create material \"%s\"; too many materials.", name));
  
  /* Note: Since a m->name is declared as a 1-element char array,
     the terminating NULL is included in sizeof(material_t) */
  m = (material_t *)malloc(sizeof(material_t)+len);
  if( m==NULL ) ERROR(("Unable to allocate material \"%s\".", name));
  m->id     = id;
  m->epsx   = epsx,   m->epsy   = epsy,   m->epsz   = epsz;
  m->mux    = mux,    m->muy    = muy,    m->muz    = muz;
  m->sigmax = sigmax, m->sigmay = sigmay, m->sigmaz = sigmaz;
  m->next   = *m_list;
  strcpy( m->name, name );
  
  *m_list = m;
  return m->id;
}

void delete_material_list( material_t **m_list ) {
  material_t *m;
  if( m_list==NULL ) return;
  while( *m_list!=NULL ) {
    m = *m_list;
    *m_list = (*m_list)->next;
    free(m);
  }
}

material_t *find_material_id( material_id id, material_t *m_list ) {
  material_t *m;
  LIST_FIND_FIRST(m,m_list,m->id==id);
  return m;
}

material_t *find_material_name( const char *name, material_t *m_list ) {
  material_t *m;
  if( name==NULL ) return NULL;
  LIST_FIND_FIRST(m,m_list,strcmp(m->name,name)==0);
  return m;
}
