#ifndef FileUtils_hxx
#define FileUtils_hxx

template<class Policy> class FileUtils_T
	: public Policy
	{
	public:

		FileUtils_T() {}
		~FileUtils_T() {}

		static int makeDirectory(const char * dirname)
			{ return Policy::makeDirectory(dirname); }

		static int getCurrentWorkingDirectory(char * dirname, size_t size)
			{ return Policy::getCurrentWorkingDirectory(dirname, size); }

	private:

	}; // class FileUtils_T

#if defined USE_MPRELAY

#if defined HOST_BUILD
#include "StandardUtilsPolicy.hxx"

typedef FileUtils_T<StandardUtilsPolicy> FileUtils;

#else
#include "P2PUtilsPolicy.hxx"

typedef FileUtils_T<P2PUtilsPolicy> FileUtils;

#endif // BUILD

#else
#include "StandardUtilsPolicy.hxx"

typedef FileUtils_T<StandardUtilsPolicy> FileUtils;

#endif

#endif // FileUtils_hxx
