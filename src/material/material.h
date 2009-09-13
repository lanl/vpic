#ifndef _material_h_
#define _material_h_

#include "../util/util.h"

enum {
  max_material = 32768 // Valid materials are numbered 0...32767
};

typedef int16_t material_id;

typedef struct material {
  char * name;                  // Name of the material
  float epsx, epsy, epsz;       // Relative permittivity along x,y,z axes
  float mux, muy, muz;          // Relative permeability along x,y,z axes
  float sigmax, sigmay, sigmaz; // Electrical conductivity along x,y,z axes
  float zetax,  zetay,  zetaz;  // Magnetic conductivity along x,y,z axes
  material_id id;               // Unique identifier for material
  struct material *next;        // Next material in list
} material_t;
  
BEGIN_C_DECLS

// In material.c

int
num_material( const material_t * m_list );

void
delete_material_list( material_t * m_list );

material_t *
find_material_id( material_id id,
                  material_t * m_list );

material_t *
find_material_name( const char * name,
                    material_t * m_list );

material_t *
append_material( material_t * m,
                 material_t ** m_list );

material_id
get_material_id( const material_t * m );

material_t *
material( const char * name,
          float epsx,   float epsy,   float epsz,
          float mux,    float muy,    float muz,
          float sigmax, float sigmay, float sigmaz,
          float zetax,  float zetay,  float zetaz );

END_C_DECLS

#endif // _material_h_
