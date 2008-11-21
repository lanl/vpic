#ifndef CheckSum_hxx
#define CheckSum_hxx

#include <openssl/evp.h>

struct CheckSum {
	unsigned char value[EVP_MAX_MD_SIZE];
	char strvalue[EVP_MAX_MD_SIZE*2+1];
	unsigned int length;
}; // struct CheckSum

template<typename T>
void md5CheckSum(T * buffer, size_t elements, CheckSum & sum) {
	size_t bytes = elements*sizeof(T);
	const EVP_MD * md = EVP_md5();

	EVP_MD_CTX ctx;

	// initialize context
	EVP_MD_CTX_init(&ctx);

	// initialize digest
	EVP_DigestInit_ex(&ctx, md, NULL);

	// update digest with buffer
	EVP_DigestUpdate(&ctx, reinterpret_cast<void *>(buffer), bytes);
	
	// finalize
	EVP_DigestFinal_ex(&ctx, sum.value, &sum.length);

	// free resources
	EVP_MD_CTX_cleanup(&ctx);

	char tmp[256];
	strcpy(sum.strvalue, "");
	for(size_t i(0); i<sum.length; i++) {
		sprintf(tmp, "%02x", sum.value[i]);
		strcat(sum.strvalue, tmp);
	} // for
} // checkSum

#endif // CheckSum_hxx
