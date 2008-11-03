/*
	Definition of FileIOData class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$
	vim: set ts=3 :
*/

#ifndef FileIOData_hxx
#define FileIOData_hxx

#include<cstddef>

enum FileIOMode {
	io_read,
	io_read_write,
	io_write,
	io_write_read,
	io_append,
	io_append_read
}; // FileIOMode

enum FileIOStatus {
	ok,
	fail
}; // FileIOStatus

#endif // FileIOData_hxx
