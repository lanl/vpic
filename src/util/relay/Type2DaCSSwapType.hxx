#ifndef Type2DaCSSwapType_hxx
#define Type2DaCSSwapType_hxx

#include <dacs.h>

template<typename T> struct Type2DaCSSwapType {};

template<> struct Type2DaCSSwapType<char>
	{
		inline static DACS_BYTE_SWAP_T type()
			{ return DACS_BYTE_SWAP_DISABLE; }
	}; // struct Type2DaCSSwapType

template<> struct Type2DaCSSwapType<int>
	{
		inline static DACS_BYTE_SWAP_T type()
			{ return DACS_BYTE_SWAP_WORD; }
	}; // struct Type2DaCSSwapType

template<> struct Type2DaCSSwapType<uint32_t>
	{
		inline static DACS_BYTE_SWAP_T type()
			{ return DACS_BYTE_SWAP_WORD; }
	}; // struct Type2DaCSSwapType

template<> struct Type2DaCSSwapType<int64_t>
	{
		inline static DACS_BYTE_SWAP_T type()
			{ return DACS_BYTE_SWAP_DOUBLE_WORD; }
	}; // struct Type2DaCSSwapType

template<> struct Type2DaCSSwapType<long>
	{
		inline static DACS_BYTE_SWAP_T type()
			{ return DACS_BYTE_SWAP_DOUBLE_WORD; }
	}; // struct Type2DaCSSwapType

template<> struct Type2DaCSSwapType<double>
	{
		inline static DACS_BYTE_SWAP_T type()
			{ return DACS_BYTE_SWAP_DOUBLE_WORD; }
	}; // struct Type2DaCSSwapType

template<> struct Type2DaCSSwapType<float>
	{
		inline static DACS_BYTE_SWAP_T type()
			{ return DACS_BYTE_SWAP_WORD; }
	}; // struct Type2DaCSSwapType

#endif // Type2DaCSSwapType_hxx
