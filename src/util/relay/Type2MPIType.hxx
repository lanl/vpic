#ifndef Type2MPIType_hxx
#define Type2MPIType_hxx

#include <mpi.h>

template<typename T> struct Type2MPIType {};

template<> struct Type2MPIType<char>
	{
		inline static MPI_Datatype type() { return MPI_BYTE; }
	}; // struct Type2MPIType

template<> struct Type2MPIType<uint32_t>
	{
		inline static MPI_Datatype type() { return MPI_INT; }
	}; // struct Type2MPIType

template<> struct Type2MPIType<int>
	{
		inline static MPI_Datatype type() { return MPI_INT; }
	}; // struct Type2MPIType

template<> struct Type2MPIType<int64_t>
	{
		inline static MPI_Datatype type() { return MPI_LONG_LONG; }
	}; // struct Type2MPIType

/*
template<> struct Type2MPIType<long long>
	{
		inline static MPI_Datatype type() { return MPI_LONG_LONG; }
	}; // struct Type2MPIType
*/

template<> struct Type2MPIType<double>
	{
		inline static MPI_Datatype type() { return MPI_DOUBLE; }
	}; // struct Type2MPIType

template<> struct Type2MPIType<float>
	{
		inline static MPI_Datatype type() { return MPI_FLOAT; }
	}; // struct Type2MPIType

#endif // Type2MPIType_hxx
