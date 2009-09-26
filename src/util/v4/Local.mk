$(call add-hdrs,v4.h v4_altivec.hxx v4_portable.hxx v4_sse.hxx)
$(call make-unit-test,test_v4,test_v4,vpicutil)
