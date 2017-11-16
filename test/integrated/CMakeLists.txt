# function to compile, run, and compare tests
function(vpic_add_test name)
  add_test(NAME ${name}
           COMMAND ${CMAKE_COMMAND}
	   -DNAME=${name}
	   -DARGS=${ARGS}
	   -DSRCDIR=${PROJECT_SOURCE_DIR}
	   -DBINDIR=${PROJECT_BINARY_DIR}
           -DMPIEXEC=${MPIEXEC}
           -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
           -DMPI_NUM_RANKS=${MPI_NUM_RANKS}
	   -DCMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}
	   -P ${PROJECT_SOURCE_DIR}/config/integrated-tests.cmake)

endfunction(vpic_add_test name)

# add the tests
set(MPI_NUM_RANKS 1)
set(ARGS "1 1")

foreach(test accel cyclo inbndj interpe outbndj pcomm)
  add_executable(${test} ${CMAKE_SOURCE_DIR}/deck/main.cc ${CMAKE_SOURCE_DIR}/deck/wrapper.cc)
  target_link_libraries(${test} vpic)
  target_compile_definitions(${test} PRIVATE INPUT_DECK=${test}.deck)
  target_include_directories(${test} PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()
vpic_add_test(cyclo)
vpic_add_test(inbndj)
vpic_add_test(interpe)
vpic_add_test(outbndj)

set(MPI_NUM_RANKS 8)

vpic_add_test(pcomm)
