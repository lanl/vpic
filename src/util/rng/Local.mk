$(call add-hdrs,rng.h)
$(call add-objs,rng rng_pool drandn_table frandn_table,vpicutil)
$(call make-unit-test,test_rng,test_rng,vpicutil)
