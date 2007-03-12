#include <spu_mfcio.h>
#include <stdio.h>
#include <test_args.h>

int main(unsigned long long speid, unsigned long long params) {
	/*
	TestArgs args;
	unsigned int tag = 0;

	spu_mfcdma32(&args, params, sizeof(TestArgs), tag, MFC_GET_CMD);
	spu_writech(MFC_WrTagMask, 1 << tag);
	(void)spu_mfcstat(2);

	printf("hello from spe %llx with id %d and value %lf\n",
		speid, args.id, args.value);

	*/
	return 0;
} // main
