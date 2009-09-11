#ifndef _material_h_
#define _material_h_

#include "../util/util.h"

/*****************************************************************************
 * Materials module
 *
 * Material properties are stored in two separate structures. A material_t is
 * a list node which contains the user specified physical material properties.
 * Every material_t is paired with a material_coefficient_t.
 * material_coefficient_t are stored in an array. The array contains the
 * numerical coefficients needed for that material in the FDTD module finite
 * difference equations. The bulk eps, mu and sigma specified in a material_t
 * are frequency indepedent. Thus, if 1 < 1/sqrt(eps*mu), the material will be
 * able to propagate signals with a superluminal group velocity. This is
 * non-physical and gives the simulation a more stringent Courant condition.
 * Typically, materials for which 1 < 1/sqrt(eps*mu) are necessarily
 * dispersive (frequency dependent) such that the group velocity is
 * subluminal. Frequency dependent materials can be modeled by adding time
 * derivative of the polarization fields to the j source term. In a similar
 * vein, conductive materials with negative dielectric constants are unstable.
 * In short, do not set eps<1 or mu<1 unless you truly understand what you are
 * doing.
 *
 * When setting the material IDs on the mesh, the material IDs should be set
 * in the ghost cells too. Further, these IDs should be consistent with the
 * neighboring domains (if any)!
 *****************************************************************************/

enum {
  max_material = 32768 // Valid materials are numbered 0...32767
};

typedef int16_t material_id;

typedef struct material {
  material_id id;               // Unique identifier for material
  float epsx, epsy, epsz;       // Relative permittivity along x,y,z axes
  float mux, muy, muz;          // Relative permeability along x,y,z axes
  float sigmax, sigmay, sigmaz; // Electrical conductivity along x,y,z axes
  float zetax,  zetay,  zetaz;  // Magnetic conductivity along x,y,z axes
  struct material *next;        // Next material in list
  char name[1];                 // Name of the material (resized at alloc)
} material_t;
  
BEGIN_C_DECLS

// In material.c

int
num_material( const material_t * m_list );

material_t *
new_material( const char * name,
              float epsx,   float epsy,   float epsz,
              float mux,    float muy,    float muz,
              float sigmax, float sigmay, float sigmaz,
              float zetax,  float zetay,  float zetaz,
              material_t ** m_list );

void
delete_material_list( material_t * m_list );

material_t *
find_material_id( material_id id,
                  material_t * m_list );

material_t *
find_material_name( const char * name,
                    material_t * m_list );

END_C_DECLS

#endif // _material_h_
