if (NO_EXPLICIT_VECTOR)
    # add the tests
    set(MPI_NUM_RANKS 1)

    set(TESTS "particle_IDs")

    # Build
    # TODO: This method of doing the tests is really bad at rebuilding them properly
    foreach(test ${TESTS})
        MESSAGE("Build")
        build_a_vpic(${test} ${CMAKE_CURRENT_SOURCE_DIR}/${test}.deck)
    endforeach()

    # Add test
    foreach(test ${TESTS})
        add_test(${test} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ./${test} ${MPIEXEC_POSTFLAGS} ${ARGS})
    endforeach()
endif(NO_EXPLICIT_VECTOR)
