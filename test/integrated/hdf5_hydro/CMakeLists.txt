# add the tests
set(MPI_NUM_RANKS 2)

set(TEST "hdf5_hydro")
find_package(Python3 COMPONENTS Interpreter)


if(USE_HDF5)
    # Build
    # TODO: This method of doing the tests is really bad at rebuilding them properly
    build_a_vpic(${TEST} ${CMAKE_CURRENT_SOURCE_DIR}/${TEST}.deck)

    # Add TEST
    add_TEST(${TEST} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPI_NUM_RANKS} ${MPIEXEC_PREFLAGS} ./${TEST} ${MPIEXEC_POSTFLAGS} ${ARGS})

    # Add a restart TEST
    list(APPEND RESTART_ARGS --restore ${CMAKE_CURRENT_BINARY_DIR}/restart0/restore.0)
    add_TEST(${TEST}_restart ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPI_NUM_RANKS} ${MPIEXEC_PREFLAGS} ./${TEST} ${MPIEXEC_POSTFLAGS} ${RESTART_ARGS})

    #and check that the output is ok
    add_TEST(${TEST}_readhdf5 python3 ${CMAKE_CURRENT_SOURCE_DIR}/readback.py ${CMAKE_CURRENT_BINARY_DIR})
endif()
