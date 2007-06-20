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

enum FileIOMode {
	io_read,
	io_write,
	io_write_append
}; // FileIOMode

enum FileIOStatus {
	ok,
	fail
}; // FileIOStatus

#endif // FileIOData_hxx
