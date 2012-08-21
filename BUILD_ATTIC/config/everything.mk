# Set any defaults not set by the including Makefile #########################

# Note:
# CC  defaults to cc
# CXX defaults to g++
# LD  defaults to ld
# AR  defaults to ar
# RM  defaults to rm -f

ifndef BUILDDIR
BUILDDIR:=default
endif

ifndef RANLIB
RANLIB:=ranlib
endif

ifndef CP
CP:=cp -f
endif

ifndef MKDIR
MKDIR:=mkdir -p
endif

ifndef SED
SED:=sed
endif

ifndef SCRUB
SCRUB:=find . -type f -name "*~" -o -name ".*.swp" -o -name "\#*" | xargs rm -fv
endif

# Explicit rules #############################################################

# make cares about these types of files
.SUFFIXES: .c .h .cxx .hxx .d .o .a

# "make all" doesn't actually make a file named "all" and similarly for others
.PHONY: all include lib bin unit-test help clean cleanall

# An unadorned "make" does "make all"
.DEFAULT_GOAL:=all

# Allows make-lib to be declared anywhere (see $$ in implicit rule below
.SECONDEXPANSION:

# Prevent make from deleting our hard won intermediate object files
.SECONDARY:

# make all should make all
all: include lib bin unit-test

# Write some diagnostic info about the configuration to the stream
help:
	# Build system configuration
	# BUILDDIR = $(BUILDDIR)
	# CC       = $(CC)
	# CXX      = $(CXX)
	# LD       = $(LD)
	# AR       = $(AR)
	# RM       = $(RM)
	# RANLIB   = $(RANLIB)
	# CP       = $(CP)
	# MKDIR    = $(MKDIR)
	# SED      = $(SED)
	# SCRUB    = $(SCRUB)
	# Explicit goals are: [include lib bin unit-test clean cleanall help]
	# "make all" does all of these but clean and help

# make clean is easy for out-of-tree builds
clean:
	$(SCRUB) && $(RM) -r objs/$(BUILDDIR)/

cleanall:
	$(SCRUB) && $(RM) -r objs

# Only do all the fancy stuff if we aren't doing a clean #####################

ifneq ($(MAKECMDGOALS),cleanall)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),help)

# Initialization #############################################################

# FIXME: This is a temporary hack to emulate the old build system from the
# end-user's point of view
EXTENSION:=$(BUILDDIR)

BUILDDIR:=objs/$(BUILDDIR)/
PREREQS:=
THISDIR:=./

# Local.mk inclusion #########################################################
 
# Finds all the Local.mk files in the source tree and includes them.  Within a
# Local.mk file is included, THISDIR is set to the path (relative to the
# Makefile dir) of the directory.  Otherwise THISDIR is the Makefile directory.

define _local-mk
$(eval THISDIR:=$(dir $(1)))
$(eval include $(1))
endef

define local-mk
$(foreach mk,$(shell find src -type f -name Local.mk),$(call _local-mk,$(mk)))
$(eval THISDIR:=./)
endef

# Local.mk helpers ###########################################################

# $(call add-objs,objs,lib)
# Adds the given objects to the given library.
# objs: Library objects.  Do not give a file extension and specify relative to
#       THISDIR.
# lib:  Library name.
define _add-objs
OBJS_$(2)+=$(foreach obj,$(1),$(BUILDDIR)$(THISDIR)$(obj).o)
PREREQS+=$(foreach obj,$(1),$(BUILDDIR)$(THISDIR)$(obj).d)
endef
add-objs = $(eval $(call _add-objs,$(1),$(2)))

# $(call add-hdrs,hdrs)
# Adds the headers the the build dir
# hdrs: Headers to add.  The file extension is necessary.  Specify relative to
#       THISDIR.
define add-hdrs
$(eval include:: $(foreach hdr,$(1),$(BUILDDIR)include/$(THISDIR)$(hdr)))
endef

# $(call make-lib,lib)
# Makes a library
# lib:  Name of the library
define make-lib
$(eval lib:: $(BUILDDIR)lib/lib$(1).a)
endef

# $(call make-bin,bin,objs,libs)
# Makes a binary
# bin:  Binary name
# objs: Binary objects.  Don not give a file extension and specify relative to
#       THISDIR
# libs: Libraries needed.
define _make-exe
$(4):: $(BUILDDIR)$(4)/$(1)
PREREQS+=$(foreach obj,$(2),$(BUILDDIR)$(THISDIR)$(obj).d)
$(BUILDDIR)$(4)/$(1): $(foreach obj,$(2),$(BUILDDIR)$(THISDIR)$(obj).o) $(foreach lib,$(3),$(BUILDDIR)lib/lib$(lib).a)
	$(MKDIR) $(BUILDDIR)$(4)
	$(LD) $(foreach obj,$(2),$(BUILDDIR)$(THISDIR)$(obj).o) $(foreach lib,$(3),$(BUILDDIR)lib/lib$(lib).a) -o $(BUILDDIR)$(4)/$(1)
endef
make-bin = $(eval $(call _make-exe,$(1),$(2),$(3),bin))

# $(call make-unit-test,bin,objs,libs)
# Makes a unit-test
# bin:  Binary name
# objs: Binary objects.  Don not give a file extension and specify relative to
#       THISDIR
# libs: Libraries needed.
make-unit-test = $(eval $(call _make-exe,$(1),$(2),$(3),unit-test))

# Implicit rules #############################################################

# How to build libraries
$(BUILDDIR)lib/lib%.a: $$(OBJS_%)
	$(MKDIR) $(dir $@) && $(RM) $@ && $(AR) cqv $@ $(OBJS_$*) && $(RANLIB) $@

# How to install headers
$(BUILDDIR)include/%:
	$(MKDIR) $(dir $@) && $(CP) $* $@

# How to build c files and their prereqs
$(BUILDDIR)%.d : %.c
	$(MKDIR) $(dir $@) && $(CC) -M $< -o $@d && $(SED) 's,^\($(notdir $(basename $@))\)\.o[ :]*,$(basename $@).o $@ : ,g' < $@d > $@ && $(RM) $@d

$(BUILDDIR)%.o : %.c
	$(MKDIR) $(dir $@) && $(CC) -c $< -o $@

# How to build c++ files and their prereqs
$(BUILDDIR)%.d : %.cxx
	$(MKDIR) $(dir $@) && $(CXX) -M $< -o $@d && $(SED) 's,^\($(notdir $(basename $@))\)\.o[ :]*,$(basename $@).o $@ : ,g' < $@d > $@ && $(RM) $@d

$(BUILDDIR)%.o : %.cxx
	$(MKDIR) $(dir $@) && $(CXX) -c $< -o $@

# Include everything #########################################################

$(eval $(local-mk))
include $(PREREQS)

# FIXME: THIS IS A TEMPORARY HACK TO EMULATE THE OLD BUILD SYSTEM FROM THE
# USER'S POINT OF VIEW

bin:: $(BUILDDIR)bin/build
$(BUILDDIR)bin/build:
	$(MKDIR) $(BUILDDIR)bin
	$(RM) $(BUILDDIR)bin/build
	echo "#!/bin/bash" > $(BUILDDIR)bin/build
	echo export INPUT_DECK_DIR='`pwd`' >> $(BUILDDIR)bin/build
	echo "pushd . >& /dev/null" >> $(BUILDDIR)bin/build
	echo cd $(PWD) >> $(BUILDDIR)bin/build
	echo $(CXX) src/main.cxx src/deck_wrapper.cxx -DINPUT_DECK='$$INPUT_DECK_DIR/$$1' -rdynamic -L$(BUILDDIR)lib -lvpic -lvpicutil -o '$$INPUT_DECK_DIR/$$1'.$(EXTENSION) >> $(BUILDDIR)bin/build
	echo "popd >& /dev/null" >> $(BUILDDIR)bin/build
	chmod u+x $(BUILDDIR)bin/build
	$(RM) build.$(EXTENSION)
	ln -s $(BUILDDIR)bin/build build.$(EXTENSION)
#$(BUILDDIR)include/src/main.cxx $(BUILDDIR)include/src/deck_wrapper.cxx

endif # Not cleanall
endif # Not clean
endif # Not help
