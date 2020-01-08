#ifndef MPWrapper_h
#define MPWrapper_h

template <class MPPolicy>
class MPWrapper_T : public MPPolicy
{
  public:
    static MPWrapper_T& instance()
    {
        static MPWrapper_T mpw;
        return mpw;
    }
    // inherited public interface
  private:
    // hide these to keep things safe
    MPWrapper_T() {}
    MPWrapper_T( const MPWrapper_T& mpw ) {}
    ~MPWrapper_T() {}
};

#ifdef USE_MPRELAY
#include "RelayPolicy.h"
typedef MPWrapper_T<RelayPolicy> MPWrapper;
#else
#include "DMPPolicy.h"
typedef MPWrapper_T<DMPPolicy> MPWrapper;
#endif

#endif // MPWrapper_h
