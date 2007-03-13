#ifndef test_args_h
#define test_args_h

#if defined __PPU__ || defined __SPU__

typedef struct
	{
		unsigned int id __attribute__ ((aligned(16)));
		double value __attribute__ ((aligned(16)));
	} TestArgs; // struct TestArgs

#endif // __PPU__ || __SPU__

#endif // test_args_h
