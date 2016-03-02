# Makefile for building final exe and cleaning directory

DECK = asymm4sp
FOURZEROSEVEN = yes
OPT = no
COMPILER = intel

# set your src path(s) and plat(form)

# defaults
SRC_PATH = /home/wohlbier/VPIC/head
PLAT = Linux
ifeq (x${FOURZEROSEVEN}, xyes)
  SRC_PATH = /home/wohlbier/VPIC/407/vpic-407-src
  PLAT = op
endif

# choose optimized or not, default is debug
FLAVOR = debug
ifeq (x${OPT}, yes)
  FLAVOR = release_sse
endif

# default to head wrapper
WRAPPER = bin/vpic
ifeq (x${FOURZEROSEVEN}, xyes)
  WRAPPER = build.${PLAT}
endif

# point to compiler wrapper
VPIC = ${SRC_PATH}/build_${COMPILER}_${FLAVOR}/${WRAPPER}

# compile final exe
all:
	$(VPIC) $(DECK)

# clean everything
clean:
	-$(RM) $(DECK).$(PLAT)
	-$(RM) -rf fields
	-$(RM) -rf global.vpc
	-$(RM) -rf hydro
	-$(RM) -rf info*
	-$(RM) -rf injectors
	-$(RM) -rf particle
	-$(RM) -rf restart*
	-$(RM) -rf rundata
	-$(RM) *~

config:
	@echo "DECK: ${DECK}"
	@echo "FOURZEROSEVEN: ${FOURZEROSEVEN}"
	@echo "OPT: ${OPT}"
	@echo "COMPILER: ${COMPILER}"
	@echo "SRC_PATH: ${SRC_PATH}"
	@echo "PLAT: ${PLAT}"
	@echo "FLAVOR: ${FLAVOR}"
	@echo "WRAPPER: ${WRAPPER}"
	@echo "VPIC: ${VPIC}"
