#ifndef vpic_unit_deck_h
#define vpic_unit_deck_h

#include "src/vpic/vpic.h"

void vpic_simulation::user_initialization( int num_cmdline_arguments,
                                           char **cmdline_arguments )
{
}

void vpic_simulation::user_diagnostics() {}

void vpic_simulation::user_particle_injection() {}

void vpic_simulation::user_current_injection() {}

void vpic_simulation::user_field_injection() {}

void vpic_simulation::user_particle_collisions() {}

#endif // vpic_unit_deck_h
