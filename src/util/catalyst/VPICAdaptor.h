#ifndef __VPICAdaptor_h
#define __VPICAdaptor_h

//#include "vpic.h"

class vpic_simulation;
class DumpParameters;
#include <string>
#include <vector>

void coprocessorinitialize (std::vector<std::string>& pythonScripts);

void coprocessorProcess (long long timestep, double time,
                         vpic_simulation* sim, int topology[3],
                         std::vector<DumpParameters *>& dumpParams);

void coprocessorfinalize ();

#endif /* __VPICAdaptor_h */
