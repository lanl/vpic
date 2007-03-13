#ifndef test_args_h
#define test_args_h

#if defined __PPU__ || defined __SPU__

/* simple example data structure for passing some information to the SPEs */
typedef struct {
	unsigned int id __attribute__ ((aligned(16)));
	double value __attribute__ ((aligned(16)));
} TestArgs;

/* simple message facility */
typedef enum {
	msg_end = 0,
	msg_advance = 1
} TestMessage;

#endif // __PPU__ || __SPU__

#endif // test_args_h
