# add the tests
set(MPI_NUM_RANKS 1)
set(ARGS "")

list(APPEND DEFUALT_ARG_TESTS simple dump)

list(APPEND RESTART_DECK dump) # Reuse existing deck and start half way
list(APPEND RESTART_BINARY restore)

foreach(test ${DEFUALT_ARG_TESTS})
    build_a_vpic(${test} ${CMAKE_CURRENT_SOURCE_DIR}/${test}.deck)
endforeach()

foreach(test ${DEFUALT_ARG_TESTS})
    add_test(${test} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${test} ${MPIEXEC_POSTFLAGS} ${ARGS})
endforeach()

# Test Restart (restore) functionality

list(APPEND CHECKPOINT_FILE "${CMAKE_CURRENT_BINARY_DIR}/checkpt.1")
list(APPEND RESTART_ARGS --restore ${CHECKPOINT_FILE})

build_a_vpic(${RESTART_BINARY} ${CMAKE_CURRENT_SOURCE_DIR}/${RESTART_DECK}.deck)
add_test(${RESTART_BINARY} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
    ${RESTART_BINARY} ${MPIEXEC_POSTFLAGS} ${RESTART_ARGS})
