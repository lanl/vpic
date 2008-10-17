/*
	Definition of StandardIOPolicy class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$
	vim: set ts=3 :
*/

#ifndef StandardIOPolicy_hxx
#define StandardIOPolicy_hxx

#include "FileIOData.hxx"
#include <string>
#include <cstdarg>
#include <cstdio>

/*!
	\class StandardIOPolicy StandardIOPolicy.h
	\brief  provides...
*/
class StandardIOPolicy
	{
	public:

		//! Constructor
		StandardIOPolicy() : handle_(NULL) {}

		//! Destructor
		~StandardIOPolicy() {}

		// open/close methods
		FileIOStatus open(const char * filename, FileIOMode mode);
		void close();
		
		int size();

		// ascii methods
		void scan(const char * format, ...);
		void print(const char * format, ...);

		// binary methods
		template<typename T> void read(T * data, size_t elements);
		template<typename T> void write(const T * data, size_t elements);

		void seek(long offset, int whence);

	private:

		FILE * handle_;

	}; // class StandardIOPolicy

inline FileIOStatus
StandardIOPolicy::open(const char * filename, FileIOMode mode)
	{
		switch(mode) {
			
			case io_read:
				handle_ = fopen(filename, "r");
				break;

			case io_read_write:
				handle_ = fopen(filename, "r+");
				break;

			case io_write:
				handle_ = fopen(filename, "w");
				break;

			case io_write_read:
				handle_ = fopen(filename, "w+");
				break;

			case io_append:
				handle_ = fopen(filename, "a");
				break;

			case io_append_read:
				handle_ = fopen(filename, "a+");
				break;

		} // switch

		if(!handle_) {
			return fail;
		} // if

		return ok;
	} // StandardIOPolicy::StandardIOPolicy

inline void StandardIOPolicy::close()
	{
		fclose(handle_);
	} // StandardIOPolicy::~StandardIOPolicy

inline int StandardIOPolicy::size()
	{
		long current = ftell(handle_);
		fseek(handle_, 0L, SEEK_END);
		int size = ftell(handle_);
		fseek(handle_, current, SEEK_SET);
		return size;
	} // StandardIOPolicy::size

inline void StandardIOPolicy::scan(const char * format, ...)
	{
		va_list ab;

		// initialize the varg list
		va_start (ab, format);

		// print to file
		vfscanf(handle_, format, ab);

		// end list
		va_end(ab);
	} // StandardIOPolicy::scan

inline void StandardIOPolicy::print(const char * format, ...)
	{
		va_list ab;

		// initialize the varg list
		va_start (ab, format);

		// print to file
		vfprintf(handle_, format, ab);

		// end list
		va_end(ab);
	} // StandardIOPolicy::print

template<typename T>
inline void StandardIOPolicy::read(T * data, size_t elements)
	{
		fread(reinterpret_cast<void *>(data), sizeof(T), elements, handle_);
	} // StandardIOPolicy::read

template<typename T>
inline void StandardIOPolicy::write(const T * data, size_t elements)
	{
		fwrite(reinterpret_cast<void *>(const_cast<T *>(data)),
			sizeof(T), elements, handle_);
	} // StandardIOPolicy::write

inline void StandardIOPolicy::seek(long offset, int whence)
	{
		fseek(handle_, offset, whence);
	} // StandardIOPolicy::seek

#endif // StandardIOPolicy_hxx
