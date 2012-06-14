#define IN_checkpt
#include <checkpt_private.h>
#include <CheckPtFileIO.hxx>

checkpt_t *
checkpt_open_rdonly( const char * name ) {
	return CheckPtFileIO::checkpt_open_rdonly(name);
}

checkpt_t *
checkpt_open_wronly( const char * name ) {
	return CheckPtFileIO::checkpt_open_wronly(name);
}

void
checkpt_close( checkpt_t * checkpt ) {
	return CheckPtFileIO::checkpt_close(checkpt);
}

void
checkpt_read( checkpt_t * checkpt,
              void * data,
              size_t sz ) {
	return CheckPtFileIO::checkpt_read(checkpt, data, sz);
}

void
checkpt_write( checkpt_t * checkpt,
               const void * data,
               size_t sz ) {
	return CheckPtFileIO::checkpt_write(checkpt, data, sz);
}
