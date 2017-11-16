# add the tests
set(MPI_NUM_RANKS 1)
set(ARGS "1 1")

foreach(test accel cyclo inbndj interpe outbndj pcomm)
  add_executable(${test} ${CMAKE_SOURCE_DIR}/deck/main.cc ${CMAKE_SOURCE_DIR}/deck/wrapper.cc)
  target_link_libraries(${test} vpic)
  target_compile_definitions(${test} PRIVATE INPUT_DECK=${test}.deck)
  target_include_directories(${test} PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()
foreach(test accel cyclo inbndj interpe outbndj)
  add_test(${test} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${test} ${MPIEXEC_POSTFLAGS} ${ARGS})
endforeach()

add_test(pcomm ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 8 ${MPIEXEC_PREFLAGS} pcomm ${MPIEXEC_POSTFLAGS} ${ARGS})
