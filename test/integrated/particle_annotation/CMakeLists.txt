# add the tests
set(MPI_NUM_RANKS 2)

set(TEST "particle_annotation")

if(PARTICLE_ANNOTATION)
    # Build
    # TODO: This method of doing the tests is really bad at rebuilding them properly
    build_a_vpic(${TEST} ${CMAKE_CURRENT_SOURCE_DIR}/${TEST}.deck)

    # Add TEST
    add_TEST(${TEST} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPI_NUM_RANKS} ${MPIEXEC_PREFLAGS} ./${TEST} ${MPIEXEC_POSTFLAGS} ${ARGS})

    # Add a restart TEST
    list(APPEND RESTART_ARGS --restore ${CMAKE_CURRENT_BINARY_DIR}/restart0/restore.0)
    add_TEST(${TEST}_restart ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPI_NUM_RANKS} ${MPIEXEC_PREFLAGS} ./${TEST} ${MPIEXEC_POSTFLAGS} ${RESTART_ARGS})
endif()
