$(call make-lib,vpicutil)
$(call add-hdrs,util_base.h util.h)
$(call add-objs,util_base boot,vpicutil)
# FIXME: EXPOSE THESE TO INPUT DECKS?
$(call add-hdrs,BitField.hxx CheckSum.hxx SystemRAM.hxx)
