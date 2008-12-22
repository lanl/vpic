#ifndef CheckSum_hxx
#define CheckSum_hxx

#if defined(ENABLE_OPENSSL)

#include <openssl/evp.h>

struct CheckSum {
	unsigned char value[EVP_MAX_MD_SIZE];
	char strvalue[EVP_MAX_MD_SIZE*2+1];
	unsigned int length;
}; // struct CheckSum

template<typename T>
void checkSumBuffer(T * buffer, size_t elements, CheckSum & sum,
	const char * digest = "md5") {
	size_t bytes = elements*sizeof(T);

	EVP_MD_CTX ctx;

	// add all digests to table
	OpenSSL_add_all_digests();

	// initialize context
	EVP_MD_CTX_init(&ctx);

	// get digest
	const EVP_MD * md = EVP_get_digestbyname(digest);
	if(!md) {
		ERROR(("Invalid digest!"));
	} // if

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
	MESSAGE(("SUM LENGTH: %d", sum.length));
	for(size_t i(0); i<sum.length; i++) {
		sprintf(tmp, "%02x", sum.value[i]);
		strcat(sum.strvalue, tmp);
	} // for
} // checkSumBuffer

#endif // ENABLE_OPENSSL

#endif // CheckSum_hxx
