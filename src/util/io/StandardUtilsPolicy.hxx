#ifndef StandardUtilsPolicy_hxx
#define StandardUtilsPolicy_hxx

#include<sys/stat.h>

class StandardUtilsPolicy
	{
	public:

		StandardUtilsPolicy() {}
		~StandardUtilsPolicy() {}

		static int makeDirectory(const char * dirname) {
			return mkdir(dirname, S_IRWXU);	
		} // mkdir

	private:

	}; // class StandardUtilsPolicy

#endif // StandardUtilsPolicy_hxx
