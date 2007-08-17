################################################################################
## VPIC library Makefile.am
################################################################################

lib_LTLIBRARIES = libvpic.la

INCLUDES = -I${top_srcdir}/src/include \
   @EXTRA_CPPFLAGS@

libvpic_includedir = ${includedir}

libvpic_include_HEADERS = \
	${top_srcdir}/src/boundary/boundary.h \
	${top_srcdir}/src/emitter/emitter.h \
	${top_srcdir}/src/field/field_pipelines.h \
	${top_srcdir}/src/field/field.h \
	${top_srcdir}/src/grid/grid.h \
	${top_srcdir}/src/include/util_base.h \
	${top_srcdir}/src/include/StandardIOPolicy.hxx \
	${top_srcdir}/src/include/mp.h \
	${top_srcdir}/src/include/species.h \
	${top_srcdir}/src/include/ConnectionManager.hxx \
	${top_srcdir}/src/include/AAISPolicy.hxx \
	${top_srcdir}/src/include/CMPolicyMultipleContextMPI.hxx \
	${top_srcdir}/src/include/v4c_spu.h \
	${top_srcdir}/src/include/FileIOData.hxx \
	${top_srcdir}/src/include/v4_spu.hxx \
	${top_srcdir}/src/include/mp_dmp.h \
	${top_srcdir}/src/include/mtrand.h \
	${top_srcdir}/src/include/vpic.hxx \
	${top_srcdir}/src/include/mp_handle.h \
	${top_srcdir}/src/include/MPWrapper.hxx \
	${top_srcdir}/src/include/emitter.h \
	${top_srcdir}/src/include/pipelines.h \
	${top_srcdir}/src/include/Type2MPIType.hxx \
	${top_srcdir}/src/include/v4_sse.hxx \
	${top_srcdir}/src/include/boundary.h \
	${top_srcdir}/src/include/DaCSUtils.h \
	${top_srcdir}/src/include/FileIO.hxx \
	${top_srcdir}/src/include/mp.hxx \
	${top_srcdir}/src/include/util.h \
	${top_srcdir}/src/include/P2PTag.hxx \
	${top_srcdir}/src/include/MPData.hxx \
	${top_srcdir}/src/include/mtrand_conv.h \
	${top_srcdir}/src/include/v4_portable.hxx \
	${top_srcdir}/src/include/P2PIOPolicy.hxx \
	${top_srcdir}/src/include/v4.h \
	${top_srcdir}/src/include/DMPConnection.hxx \
	${top_srcdir}/src/include/MPRelay.hxx \
	${top_srcdir}/src/include/grid.h \
	${top_srcdir}/src/include/field_pipelines.h \
	${top_srcdir}/src/include/P2PPolicyDaCS.hxx \
	${top_srcdir}/src/include/particle.h \
	${top_srcdir}/src/include/CMPolicyDaCS.hxx \
	${top_srcdir}/src/include/CMPolicyMPIDaCS.hxx \
	${top_srcdir}/src/include/mp_t.h \
	${top_srcdir}/src/include/field.h \
	${top_srcdir}/src/include/Type2DaCSSwapType.hxx \
	${top_srcdir}/src/include/P2PPolicyMPI.hxx \
	${top_srcdir}/src/include/DMPPolicy.hxx \
	${top_srcdir}/src/include/swap.h \
	${top_srcdir}/src/include/particle_pipelines.h \
	${top_srcdir}/src/include/material.h \
	${top_srcdir}/src/include/P2PConnection.hxx \
	${top_srcdir}/src/include/DMPPolicyMPI.hxx \
	${top_srcdir}/src/material/material.h \
	${top_srcdir}/src/particle/particle.h \
	${top_srcdir}/src/particle/particle_pipelines.h \
	${top_srcdir}/src/species/species.h \
	${top_srcdir}/src/util/util_base.h \
	${top_srcdir}/src/util/relay/ConnectionManager.hxx \
	${top_srcdir}/src/util/relay/CMPolicyMultipleContextMPI.hxx \
	${top_srcdir}/src/util/relay/Type2MPIType.hxx \
	${top_srcdir}/src/util/relay/p2p/DaCSUtils.h \
	${top_srcdir}/src/util/relay/p2p/P2PTag.hxx \
	${top_srcdir}/src/util/relay/p2p/P2PPolicyDaCS.hxx \
	${top_srcdir}/src/util/relay/p2p/P2PPolicyMPI.hxx \
	${top_srcdir}/src/util/relay/p2p/P2PConnection.hxx \
	${top_srcdir}/src/util/relay/MPData.hxx \
	${top_srcdir}/src/util/relay/dmp/DMPConnection.hxx \
	${top_srcdir}/src/util/relay/dmp/DMPPolicyMPI.hxx \
	${top_srcdir}/src/util/relay/MPRelay.hxx \
	${top_srcdir}/src/util/relay/CMPolicyDaCS.hxx \
	${top_srcdir}/src/util/relay/CMPolicyMPIDaCS.hxx \
	${top_srcdir}/src/util/relay/Type2DaCSSwapType.hxx \
	${top_srcdir}/src/util/pipelines/pipelines.h \
	${top_srcdir}/src/util/mp/mp.h \
	${top_srcdir}/src/util/mp/aais/AAISPolicy.hxx \
	${top_srcdir}/src/util/mp/mp_handle.h \
	${top_srcdir}/src/util/mp/MPWrapper.hxx \
	${top_srcdir}/src/util/mp/mp.hxx \
	${top_srcdir}/src/util/mp/dmp/mp_dmp.h \
	${top_srcdir}/src/util/mp/dmp/mp_t.h \
	${top_srcdir}/src/util/mp/dmp/DMPPolicy.hxx \
	${top_srcdir}/src/util/mtrand/mtrand.h \
	${top_srcdir}/src/util/mtrand/mtrand_conv.h \
	${top_srcdir}/src/util/util.h \
	${top_srcdir}/src/util/io/StandardIOPolicy.hxx \
	${top_srcdir}/src/util/io/FileIOData.hxx \
	${top_srcdir}/src/util/io/FileIO.hxx \
	${top_srcdir}/src/util/io/P2PIOPolicy.hxx \
	${top_srcdir}/src/util/v4/v4c_spu.h \
	${top_srcdir}/src/util/v4/v4_spu.hxx \
	${top_srcdir}/src/util/v4/v4_sse.hxx \
	${top_srcdir}/src/util/v4/v4_portable.hxx \
	${top_srcdir}/src/util/v4/v4.h \
	${top_srcdir}/src/util/swap.h \
	${top_srcdir}/src/vpic/vpic.hxx

libvpic_la_SOURCES = \
	${top_srcdir}/src/boundary/maxwellian_reflux.c \
	${top_srcdir}/src/boundary/absorb_tally.c \
	${top_srcdir}/src/boundary/link.c \
	${top_srcdir}/src/emitter/child-langmuir.c \
	${top_srcdir}/src/emitter/ivory.c \
	${top_srcdir}/src/emitter/ccube.c \
	${top_srcdir}/src/emitter/emitter.c \
	${top_srcdir}/src/field/compute_div_b_err.cxx \
	${top_srcdir}/src/field/compute_div_e_err.cxx \
	${top_srcdir}/src/field/distribute_voxels.c \
	${top_srcdir}/src/field/clean_div_b.cxx \
	${top_srcdir}/src/field/remote.c \
	${top_srcdir}/src/field/compute_rhob.cxx \
	${top_srcdir}/src/field/field_structors.c \
	${top_srcdir}/src/field/energy_f.cxx \
	${top_srcdir}/src/field/unload_accumulator.cxx \
	${top_srcdir}/src/field/compute_curl_b.cxx \
	${top_srcdir}/src/field/reduce_accumulators.cxx \
	${top_srcdir}/src/field/advance_b.cxx \
	${top_srcdir}/src/field/clear_accumulators.cxx \
	${top_srcdir}/src/field/clean_div_e.cxx \
	${top_srcdir}/src/field/local.c \
	${top_srcdir}/src/field/load_interpolator.cxx \
	${top_srcdir}/src/field/compute_rms_div_e_err.cxx \
	${top_srcdir}/src/field/compute_rms_div_b_err.cxx \
	${top_srcdir}/src/field/advance_e.cxx \
	${top_srcdir}/src/grid/ops.c \
	${top_srcdir}/src/grid/add_boundary.c \
	${top_srcdir}/src/grid/partition.c \
	${top_srcdir}/src/grid/grid_structors.c \
	${top_srcdir}/src/material/coefficient.c \
	${top_srcdir}/src/material/material.c \
	${top_srcdir}/src/particle/uncenter_p.cxx \
	${top_srcdir}/src/particle/boundary_p.c \
	${top_srcdir}/src/particle/particle_structors.c \
	${top_srcdir}/src/particle/advance_p.cxx \
	${top_srcdir}/src/particle/sort_p.c \
	${top_srcdir}/src/particle/move_p.c \
	${top_srcdir}/src/particle/center_p.cxx \
	${top_srcdir}/src/particle/energy_p.cxx \
	${top_srcdir}/src/particle/rho_p.c \
	${top_srcdir}/src/particle/hydro_p.c \
	${top_srcdir}/src/species/species.c \
	${top_srcdir}/src/util/util.c \
	${top_srcdir}/src/util/pipelines/pipelines_spu.c \
	${top_srcdir}/src/util/pipelines/pipelines_serial.c \
	${top_srcdir}/src/util/pipelines/pipelines_thread.c \
	${top_srcdir}/src/util/mp/dmp/mp_dmp.c \
	${top_srcdir}/src/util/mp/mp.cxx \
	${top_srcdir}/src/util/mtrand/mtrand.c \
	${top_srcdir}/src/vpic/vpic.cxx \
	${top_srcdir}/src/vpic/advance.cxx \
	${top_srcdir}/src/vpic/initialize.cxx \
	${top_srcdir}/src/vpic/misc.cxx \
	${top_srcdir}/src/vpic/dump.cxx

libvpic_la_LDFLAGS = @VERSION_INFORMATION@ \
	@EXTRA_LDFLAGS@
libvpic_la_LIBADD = @EXTRA_LIBS@
