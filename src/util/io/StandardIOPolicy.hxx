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
		
		bool isOpen() { return is_open_; }

		// return file size in bytes
		uint64_t size();

		// ascii methods
		void print(const char * format, va_list & args);

		// binary methods
		template<typename T> void read(T * data, size_t elements);
		template<typename T> void write(const T * data, size_t elements);

		void seek(uint64_t offset, int32_t whence);
		uint64_t tell();
		void rewind();

	private:

		bool is_open_;
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

		is_open_ = true;
		return ok;
	} // StandardIOPolicy::StandardIOPolicy

inline void StandardIOPolicy::close()
	{
		fclose(handle_);
		is_open_ = false;
	} // StandardIOPolicy::~StandardIOPolicy

inline uint64_t StandardIOPolicy::size()
	{
		long current = ftell(handle_);
		fseek(handle_, 0L, SEEK_END);
		int size = ftell(handle_);
		fseek(handle_, current, SEEK_SET);
		return (uint64_t)size;
	} // StandardIOPolicy::size

inline void StandardIOPolicy::print(const char * format, va_list & args)
	{
		// print to file
		vfprintf(handle_, format, args);

		// end list
		va_end(args);
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

inline void StandardIOPolicy::seek(uint64_t offset, int32_t whence)
	{
		fseek(handle_, offset, whence);
	} // StandardIOPolicy::seek

inline uint64_t StandardIOPolicy::tell()
	{
		return uint64_t(ftell(handle_));
	} // StandardIOPolicy::tell

inline void StandardIOPolicy::rewind()
	{
		StandardIOPolicy::seek(uint64_t(0), SEEK_SET);
	} // StandardIOPolicy::tell

#endif // StandardIOPolicy_hxx
