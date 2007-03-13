#if defined __SPU__

#include <spu_mfcio.h>
#include <stdio.h>
#include <test_args.h>

typedef union {
	unsigned long long ull;
	unsigned int ui[2];
} addr_t;

int main(unsigned long long speid, unsigned long long params) {
	TestArgs args __attribute__ ((aligned(16)));
	unsigned int tag = 0;
	TestMessage msg;
	addr_t addr;

	addr.ull = params;

	spu_mfcdma32(&args, addr.ui[0], sizeof(TestArgs), tag, MFC_GET_CMD);
	spu_writech(MFC_WrTagMask, 1 << tag);
	(void)spu_mfcstat(2);

	printf("hello from spe %llx with id %d and value %lf\n",
		speid, args.id, args.value);

	do {
		msg = spu_read_in_mbox();

		switch(msg) {
			case msg_advance:
				printf("spe %llx with id %d was told to advance\n",
					speid, args.id);
				break;
		} // switch

	} while(msg != msg_end);

	return 0;
} // main

#endif // __SPU__
