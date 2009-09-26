$(call add-hdrs,sf_interface.h)
# FIXME: MERGE ACCUMULATORS OBJECTS
$(call add-objs,interpolator_array accumulator_array clear_accumulators reduce_accumulators unload_accumulator hydro_array,vpic)
