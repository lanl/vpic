# add the tests
set(MPI_NUM_RANKS 1)
set(ARGS "1 1")

foreach(test accel cyclo inbndj interpe outbndj pcomm)
  build_a_vpic(${test} ${CMAKE_CURRENT_SOURCE_DIR}/${test}.deck)
endforeach()
foreach(test accel cyclo inbndj interpe outbndj)
  add_test(${test} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${test} ${MPIEXEC_POSTFLAGS} ${ARGS})
endforeach()

add_test(pcomm ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 8 ${MPIEXEC_PREFLAGS} pcomm ${MPIEXEC_POSTFLAGS} ${ARGS})
