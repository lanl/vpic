# Makefile for building final exe and cleaning directory

# name of platform to build for
#PLAT = op
PLAT = emt64

# path to installed compiler wrapper
#VPIC = /home/wohlbier/VPIC/407/vpic-407-src/build_intel_debug/build.${PLAT}
VPIC = /home/wohlbier/VPIC/407/vpic-407-src/build_intel_release_sse/build.${PLAT}

# name of input deck file
DECK = open-collisional

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

# example openmpi mpirun commands. --tpp is number of vpic pipelines. Each
# pipeline is domain replicated.

# one rank on a node, bound to all available cores.
# mpirun bind-to none  --report-bindings -n 1 ./harris.Linux  1 1 --tpp 16

# two ranks on a node with two sockets. ranks bound to sockets.
# mpirun --map-by socket --bind-to socket --report-bindings -n 2 ./harris.Linux  1 1 --tpp 16

# 4 ranks on a node with 4 numa nodes. ranks bound to numa nodes.
# mpirun --map-by numa --bind-to numa --report-bindings -n 4 ./harris.Linux  1 1 --tpp 8
