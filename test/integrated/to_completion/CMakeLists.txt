# WARNING: None of these tests test correctness, only that they don't die.

set(MPIEXEC_NUMPROC 1)

# TODO: If would be good if we could detect that the machine supports the
# desired number of threads
set(MPIEXEC_NUMPROC_PARALLEL 8)

set(ARGS "")

list(APPEND DEFUALT_ARG_TESTS
    simple # This test is a simple run which should not die
    dump) # This is a simple run which should dump restart files

list(APPEND RESTART_DECK dump) # Reuse existing deck and start half way
list(APPEND RESTART_BINARY restore)

foreach(test ${DEFUALT_ARG_TESTS})
    build_a_vpic(${test} ${CMAKE_CURRENT_SOURCE_DIR}/${test}.deck)
endforeach()

foreach(test ${DEFUALT_ARG_TESTS})
    add_test(${test} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_NUMPROC}
        ${MPIEXEC_PREFLAGS} ${test} ${MPIEXEC_POSTFLAGS} ${ARGS}) endforeach()

# Try a parallel run
set (PARALLEL_TEST parallel)
build_a_vpic(${PARALLEL_TEST} ${CMAKE_CURRENT_SOURCE_DIR}/simple.deck)
add_test(${PARALLEL_TEST} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG}
    ${MPIEXEC_NUMPROC_PARALLEL} ${MPIEXEC_PREFLAGS} ${PARALLEL_TEST}
    ${MPIEXEC_POSTFLAGS} ${ARGS})

# Try a threaded run
set (THREADED_TEST threaded)
list(APPEND THREADED_ARGS --tpp ${MPIEXEC_NUMPROC_PARALLEL})

build_a_vpic(${THREADED_TEST} ${CMAKE_CURRENT_SOURCE_DIR}/simple.deck)
add_test(${THREADED_TEST} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_NUMPROC}
    ${MPIEXEC_PREFLAGS} ${THREADED_TEST} ${MPIEXEC_POSTFLAGS} ${THREADED_ARGS})

# TODO: Do we want to try an MPI + Threaded runs

# Test Restart (restore) functionality

list(APPEND CHECKPOINT_FILE "${CMAKE_CURRENT_BINARY_DIR}/checkpt.1")
list(APPEND RESTART_ARGS --restore ${CHECKPOINT_FILE})

build_a_vpic(${RESTART_BINARY} ${CMAKE_CURRENT_SOURCE_DIR}/${RESTART_DECK}.deck)
add_test(${RESTART_BINARY} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG}
    ${MPIEXEC_NUMPROC} ${MPIEXEC_PREFLAGS} ${RESTART_BINARY}
    ${MPIEXEC_POSTFLAGS} ${RESTART_ARGS})
