#ifndef test_args_h
#define test_args_h

struct TestArgs
	{
		unsigned int id __attribute__ ((aligned(16)));
		double value __attribute__ ((aligned(16)));
	}; // class TestArgs

#endif // test_args_h
