/*------------------------------------------------------------------------------
	Wrapper class to facilitate changing the message passing interface
	without having to do crazy compiler stuff.

	The original C implementation is in the "dmp" sub-directory for
	distributed memory parallel...

	The "aais" directory will contain an Opteron <--> Cell style
	implementation.
------------------------------------------------------------------------------*/
#ifndef MPWrapper_hxx
#define MPWrapper_hxx

#include <DMPPolicy.hxx>
#include <AAISPolicy.hxx>

template<class MPPolicy> class MPWrapper_T
	: public MPPolicy
	{
	public:

		// Myer's singleton just in case this class
		// needs to maintain state information
		// for some policy implementations
		static MPWrapper_T & instance() {
			static MPWrapper_T mpw;
			return mpw;
		} // instance

		// inherited public interface
	
	private:

		// hide these to keep things safe
		MPWrapper_T() {}
		MPWrapper_T(const MPWrapper_T & mpw) {}
		~MPWrapper_T() {}

	}; // class MPWrapper_T

// compile time selection of implementation policy
#if defined USE_AAIS_MP
	typedef MPWrapper_T<AAISPolicy> MPWrapper;
#else
	typedef MPWrapper_T<DMPPolicy> MPWrapper;
#endif // MP Implementation

#endif // MPWrapper_hxx
