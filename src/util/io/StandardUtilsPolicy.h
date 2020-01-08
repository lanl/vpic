#ifndef StandardUtilsPolicy_h
#define StandardUtilsPolicy_h

#include <sys/stat.h>
#include <unistd.h>

class StandardUtilsPolicy
{
  public:
    StandardUtilsPolicy() {}
    ~StandardUtilsPolicy() {}

    static int makeDirectory( const char *dirname )
    {
        return mkdir( dirname, S_IRWXU );
    } // makeDirectory

    static int getCurrentWorkingDirectory( char *dirname, size_t size )
    {
        if ( getcwd( dirname, size ) == NULL )
        {
            return -1;
        }
        else
        {
            return 0;
        }
    } // getCurrentWorkingDirectory

  private:
}; // class StandardUtilsPolicy

#endif // StandardUtilsPolicy_h
